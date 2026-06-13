function val = ibl_validate_within_session_stitching(qc, cfg, varargin)
%IBL_VALIDATE_WITHIN_SESSION_STITCHING Split-half validation of CSD stitching.
%
%   val = ibl_validate_within_session_stitching(qc, cfg)
%
% The first half of each session is used to create artificial partial
% pairwise observations. Held-out pairs are completed by the stitcher and
% compared against the full simultaneous CSD from the second half. The full
% first-half simultaneous CSD is reported as the split-half upper bound.

assert(exist('CrossSpecMatpt', 'file') == 2, ...
    'Chronux CrossSpecMatpt not found. Add Chronux to the MATLAB path.');

p = inputParser;
p.addParameter('SessionIds', strings(0,1), @(x)isstring(x)||iscellstr(x)||ischar(x));
p.addParameter('MaxSessions', 3, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('MinAreas', 6, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('MaxAreas', 12, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('MinDuration', 600, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('HoldoutFraction', 0.25, @(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
p.addParameter('Seed', 1, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('ComputePDC', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('MakeFigure', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('OutputFile', "", @(x)isstring(x)||ischar(x));
p.parse(varargin{:});
opt = p.Results;
opt.SessionIds = string(opt.SessionIds);
if strlength(string(opt.OutputFile)) == 0
    outDir = fullfile(cfg.output_root, 'validation');
    if exist(outDir, 'dir') ~= 7
        mkdir(outDir);
    end
    opt.OutputFile = fullfile(outDir, 'within_session_stitching.mat');
end
outDir = fileparts(char(opt.OutputFile));
if ~isempty(outDir) && exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end

if nargin < 1 || isempty(qc)
    tmp = load(cfg.pca_qc_file, 'qc');
    qc = tmp.qc;
end

rng(opt.Seed);
sessions = local_select_sessions(qc, cfg, opt);
records = struct([]);

if cfg.verbose
    fprintf('[ibl_validate_within_session_stitching] validating %d sessions\n', numel(sessions));
end

for k = 1:numel(sessions)
    sessionId = sessions(k).session_id;
    if cfg.verbose
        fprintf('[within-session validation] %d/%d %s\n', k, numel(sessions), sessionId);
    end

    pcaFile = fullfile(cfg.area_pca_dir, [sessionId '.mat']);
    tmp = load(pcaFile, 'summary');
    pcaSummary = tmp.summary;

    areaNames = local_select_areas(pcaSummary, qc.largest_component, opt);
    if numel(areaNames) < opt.MinAreas
        if cfg.verbose
            fprintf('  skip: only %d selected areas\n', numel(areaNames));
        end
        continue;
    end

    halfDur = pcaSummary.sp_dur / 2;
    [Strain, f] = local_compute_window_csd(pcaSummary, cfg, areaNames, [0 halfDur]);
    [Stest, f2] = local_compute_window_csd(pcaSummary, cfg, areaNames, [halfDur pcaSummary.sp_dur]);
    assert(numel(f) == numel(f2) && max(abs(f(:)-f2(:))) < 1e-10, ...
        'Split-half frequency mismatch for %s.', sessionId);

    nAreas = numel(areaNames);
    holdoutMask = local_holdout_mask(nAreas, opt.HoldoutFraction);
    [Sblocks, recset] = local_pair_blocks(Strain, holdoutMask);

    stitchParams = cfg.stitch;
    stitchParams.verbose = false;
    [~, Sstitch, SregStitch, f, meta] = stitch_spectra_blocks(Sblocks, recset, f, stitchParams);

    rec = struct();
    rec.session_id = sessionId;
    rec.area_names = areaNames;
    rec.freqs = f;
    rec.holdout_mask = holdoutMask;
    rec.observed_pair_count = meta.count_mat;
    rec.metrics = local_csd_metrics(Sstitch, Strain, Stest, holdoutMask);

    if opt.ComputePDC
        pdcParams = stitchParams;
        pdcParams.regularizer = 'eigfloor';
        pdcParams.lambda = 0;
        [~, SregTrain] = glasso_precision_estimate(Strain, 0, pdcParams);
        [~, SregTest] = glasso_precision_estimate(Stest, 0, pdcParams);
        pdcStitch = nonparam_pdc_H(SregStitch, f, ...
            'maxiter', cfg.chronux.maxiter, 'tol', cfg.chronux.tol);
        pdcTrain = nonparam_pdc_H(SregTrain, f, ...
            'maxiter', cfg.chronux.maxiter, 'tol', cfg.chronux.tol);
        pdcTest = nonparam_pdc_H(SregTest, f, ...
            'maxiter', cfg.chronux.maxiter, 'tol', cfg.chronux.tol);
        rec.infoflow_stitch = ibl_pdc_to_infoflow(pdcStitch);
        rec.infoflow_train_upper = ibl_pdc_to_infoflow(pdcTrain);
        rec.infoflow_test = ibl_pdc_to_infoflow(pdcTest);
        rec.infoflow_metrics = local_infoflow_metrics( ...
            rec.infoflow_stitch, rec.infoflow_train_upper, rec.infoflow_test, holdoutMask);
    end

    if opt.MakeFigure
        local_plot_within_record(rec);
    end

    records = [records; rec]; %#ok<AGROW>
    if cfg.verbose
        fprintf('  heldout CSD abs corr stitch=%.3f upper=%.3f  relRMSE stitch=%.3f upper=%.3f\n', ...
            rec.metrics.heldout_abs_corr_stitch, rec.metrics.heldout_abs_corr_upper, ...
            rec.metrics.heldout_rel_rmse_stitch, rec.metrics.heldout_rel_rmse_upper);
        if isfield(rec, 'infoflow_metrics')
            fprintf('  heldout infoflow corr stitch=%.3f upper=%.3f\n', ...
                rec.infoflow_metrics.heldout_corr_stitch, rec.infoflow_metrics.heldout_corr_upper);
        end
    end
end

val.kind = "within_session_split_half";
val.records = records;
val.options = opt;
val.created_at = string(datetime('now'));
save(opt.OutputFile, 'val', '-v7.3');
if cfg.verbose
    fprintf('[ibl_validate_within_session_stitching] saved %s\n', opt.OutputFile);
end
end


function sessions = local_select_sessions(qc, cfg, opt)
sessions = qc.qualifying_sessions;
if ~isempty(opt.SessionIds)
    keep = ismember(string({sessions.session_id}), opt.SessionIds);
    sessions = sessions(keep);
end

selected = struct([]);
for k = 1:numel(sessions)
    pcaFile = fullfile(cfg.area_pca_dir, [sessions(k).session_id '.mat']);
    if exist(pcaFile, 'file') ~= 2
        continue;
    end
    tmp = load(pcaFile, 'summary');
    areaNames = local_select_areas(tmp.summary, qc.largest_component, opt);
    if tmp.summary.sp_dur >= opt.MinDuration && numel(areaNames) >= opt.MinAreas
        selected = [selected; sessions(k)]; %#ok<AGROW>
    end
    if numel(selected) >= opt.MaxSessions
        break;
    end
end
sessions = selected;
end


function areaNames = local_select_areas(pcaSummary, allowedAreas, opt)
scores = pcaSummary.mean_pc1_explained(:);
keep = pcaSummary.pass_qc_area(:) & ismember(string(pcaSummary.area_names(:)), allowedAreas) & ...
    isfinite(scores);
areaNames = string(pcaSummary.area_names(keep));
scores = scores(keep);
[~, ord] = sort(scores, 'descend');
ord = ord(1:min(numel(ord), opt.MaxAreas));
areaNames = sort(areaNames(ord));
end


function [Sarea, f] = local_compute_window_csd(pcaSummary, cfg, areaNames, win)
regionTable = [];
if isfield(cfg, 'region_table')
    regionTable = cfg.region_table;
end
if isempty(regionTable) && isfield(pcaSummary, 'region_table')
    regionTable = pcaSummary.region_table;
end
if isempty(regionTable)
    tmp = load(cfg.scan_file, 'scan');
    regionTable = tmp.scan.region_table;
end

areaStruct = ibl_load_area_spikes(pcaSummary.session_root, regionTable, cfg);
[allNames, ord] = sort(string({areaStruct.area}));
areaStruct = areaStruct(ord);
[tf, idxStruct] = ismember(areaNames, allNames);
assert(all(tf), 'Could not reload all selected areas for %s.', pcaSummary.session_id);

[~, idxPca] = ismember(areaNames, string(pcaSummary.area_names(:)));
fixedLoadings = pcaSummary.fixed_loadings(idxPca);

Tmax = win(2) - win(1);
nAreas = numel(areaNames);
Sarea = [];
f = [];

for a = 1:nAreas
    data = local_times_to_struct(areaStruct(idxStruct(a)).spike_times, win);
    [Saa, fa] = local_point_spectrum(data, Tmax, cfg.chronux);
    [Saa, fa] = ibl_reduce_frequency_grid(Saa, fa, cfg.target_n_freqs);
    if isempty(f)
        f = fa(:);
        Sarea = zeros(nAreas, nAreas, numel(f));
    else
        assert(numel(f) == numel(fa) && max(abs(f(:)-fa(:))) < 1e-10, ...
            'Frequency mismatch while computing validation auto spectra.');
    end
    va = fixedLoadings{a};
    Sarea(a,a,:) = ibl_project_cross_spectrum_fixed(Saa, va, va);
end

for a = 1:nAreas
    for b = (a+1):nAreas
        data = [local_times_to_struct(areaStruct(idxStruct(a)).spike_times, win); ...
                local_times_to_struct(areaStruct(idxStruct(b)).spike_times, win)];
        [SabFull, fab] = local_point_spectrum(data, Tmax, cfg.chronux);
        [SabFull, fab] = ibl_reduce_frequency_grid(SabFull, fab, cfg.target_n_freqs);
        assert(numel(f) == numel(fab) && max(abs(f(:)-fab(:))) < 1e-10, ...
            'Frequency mismatch while computing validation cross spectra.');
        na = areaStruct(idxStruct(a)).n_neurons;
        nb = areaStruct(idxStruct(b)).n_neurons;
        Sab = SabFull(1:na, na+1:na+nb, :);
        vab = ibl_project_cross_spectrum_fixed(Sab, fixedLoadings{a}, fixedLoadings{b});
        Sarea(a,b,:) = vab;
        Sarea(b,a,:) = conj(vab);
    end
end
end


function data = local_times_to_struct(timesCell, win)
n = numel(timesCell);
data = struct('times', cell(n,1));
for i = 1:n
    t = timesCell{i}(:);
    t = t(t >= win(1) & t < win(2)) - win(1);
    data(i).times = t;
end
end


function [S, f] = local_point_spectrum(data, Tmax, chronux)
[pxy, ~, ~, ~, ~, f] = CrossSpecMatpt(data, chronux.win, Tmax, chronux);
S = permute(pxy, [3, 2, 1]);
for i = 1:size(S,3)
    S(:,:,i) = (S(:,:,i) + S(:,:,i)') / 2;
end
S(~isfinite(S)) = 0;
end


function mask = local_holdout_mask(nAreas, frac)
pairs = nchoosek(1:nAreas, 2);
nHold = max(1, round(frac * size(pairs, 1)));
ord = randperm(size(pairs, 1), nHold);
mask = false(nAreas, nAreas);
for k = 1:numel(ord)
    i = pairs(ord(k), 1);
    j = pairs(ord(k), 2);
    mask(i,j) = true;
    mask(j,i) = true;
end
end


function [Sblocks, recset] = local_pair_blocks(S, holdoutMask)
nAreas = size(S, 1);
Sblocks = {};
recset = {};
for i = 1:nAreas
    Sblocks{end+1} = S(i,i,:); %#ok<AGROW>
    recset{end+1} = i; %#ok<AGROW>
end
for i = 1:nAreas
    for j = (i+1):nAreas
        if holdoutMask(i,j)
            continue;
        end
        Sblocks{end+1} = S([i j], [i j], :); %#ok<AGROW>
        recset{end+1} = [i j]; %#ok<AGROW>
    end
end
end


function metrics = local_csd_metrics(Sstitch, Strain, Stest, holdoutMask)
metrics.heldout_abs_corr_stitch = local_abs_corr(Sstitch, Stest, holdoutMask);
metrics.heldout_abs_corr_upper = local_abs_corr(Strain, Stest, holdoutMask);
metrics.heldout_rel_rmse_stitch = local_rel_rmse(Sstitch, Stest, holdoutMask);
metrics.heldout_rel_rmse_upper = local_rel_rmse(Strain, Stest, holdoutMask);
obsMask = ~holdoutMask & ~eye(size(holdoutMask));
metrics.observed_abs_corr_stitch = local_abs_corr(Sstitch, Stest, obsMask);
metrics.observed_abs_corr_upper = local_abs_corr(Strain, Stest, obsMask);
metrics.n_heldout_pairs = nnz(triu(holdoutMask, 1));
metrics.n_observed_pairs = nnz(triu(obsMask, 1));
end


function metrics = local_infoflow_metrics(Fstitch, Ftrain, Ftest, holdoutMask)
dirMask = holdoutMask & ~eye(size(holdoutMask));
obsMask = ~holdoutMask & ~eye(size(holdoutMask));
metrics.heldout_corr_stitch = local_vec_corr(Fstitch(dirMask), Ftest(dirMask));
metrics.heldout_corr_upper = local_vec_corr(Ftrain(dirMask), Ftest(dirMask));
metrics.observed_corr_stitch = local_vec_corr(Fstitch(obsMask), Ftest(obsMask));
metrics.observed_corr_upper = local_vec_corr(Ftrain(obsMask), Ftest(obsMask));
metrics.heldout_rel_rmse_stitch = norm(Fstitch(dirMask) - Ftest(dirMask)) / max(norm(Ftest(dirMask)), eps);
metrics.heldout_rel_rmse_upper = norm(Ftrain(dirMask) - Ftest(dirMask)) / max(norm(Ftest(dirMask)), eps);
end


function r = local_abs_corr(A, B, pairMask)
idx = find(triu(pairMask, 1));
if isempty(idx)
    r = NaN;
    return;
end
nf = size(A, 3);
va = zeros(numel(idx) * nf, 1);
vb = zeros(numel(idx) * nf, 1);
pos = 1;
for k = 1:numel(idx)
    [i, j] = ind2sub(size(pairMask), idx(k));
    a = squeeze(abs(A(i,j,:)));
    b = squeeze(abs(B(i,j,:)));
    va(pos:pos+nf-1) = a(:);
    vb(pos:pos+nf-1) = b(:);
    pos = pos + nf;
end
r = local_vec_corr(va, vb);
end


function e = local_rel_rmse(A, B, pairMask)
idx = find(triu(pairMask, 1));
if isempty(idx)
    e = NaN;
    return;
end
va = [];
vb = [];
for k = 1:numel(idx)
    [i, j] = ind2sub(size(pairMask), idx(k));
    va = [va; squeeze(A(i,j,:))]; %#ok<AGROW>
    vb = [vb; squeeze(B(i,j,:))]; %#ok<AGROW>
end
e = norm(va - vb) / max(norm(vb), eps);
end


function r = local_vec_corr(a, b)
a = real(a(:));
b = real(b(:));
good = isfinite(a) & isfinite(b);
if nnz(good) < 3 || std(a(good)) == 0 || std(b(good)) == 0
    r = NaN;
else
    r = corr(a(good), b(good), 'type', 'Pearson');
end
end


function local_plot_within_record(rec)
figure('Color', 'w', 'Position', [100 100 1180 430]);
tiledlayout(1,3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
imagesc(rec.holdout_mask);
axis image;
title(sprintf('%s held-out pairs', rec.session_id), 'Interpreter', 'none');
set(gca, 'XTick', 1:numel(rec.area_names), 'XTickLabel', rec.area_names, ...
    'YTick', 1:numel(rec.area_names), 'YTickLabel', rec.area_names, ...
    'XTickLabelRotation', 90, 'TickLabelInterpreter', 'none');

nexttile;
bar([rec.metrics.heldout_abs_corr_stitch, rec.metrics.heldout_abs_corr_upper; ...
     rec.metrics.heldout_rel_rmse_stitch, rec.metrics.heldout_rel_rmse_upper]);
set(gca, 'XTickLabel', {'abs corr', 'rel RMSE'});
legend({'stitched', 'upper'}, 'Location', 'best');
title('Held-out CSD recovery');

nexttile;
if isfield(rec, 'infoflow_metrics')
    bar([rec.infoflow_metrics.heldout_corr_stitch, rec.infoflow_metrics.heldout_corr_upper; ...
         rec.infoflow_metrics.heldout_rel_rmse_stitch, rec.infoflow_metrics.heldout_rel_rmse_upper]);
    set(gca, 'XTickLabel', {'corr', 'rel RMSE'});
    legend({'stitched', 'upper'}, 'Location', 'best');
    title('Held-out infoflow recovery');
else
    axis off;
    text(0.1, 0.5, 'PDC/infoflow skipped');
end
end
