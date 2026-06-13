function val = ibl_validate_within_session_subset_stitching(qc, cfg, varargin)
%IBL_VALIDATE_WITHIN_SESSION_SUBSET_STITCHING Synthetic partial-observation validation.
%
%   val = ibl_validate_within_session_subset_stitching(qc, cfg)
%
% For each selected session, this function creates multiple time-separated
% partial observations. Each observation sees only a subset of areas, subsets
% overlap, and the union covers all selected areas. Area PC loadings are
% computed independently inside each observation window, matching the real
% across-session use case. The stitched CSD is compared against a full
% simultaneous-observation control computed from all selected areas in the
% same session.

assert(exist('CrossSpecMatpt', 'file') == 2, ...
    'Chronux CrossSpecMatpt not found. Add Chronux to the MATLAB path.');

p = inputParser;
p.addParameter('SessionIds', strings(0,1), @(x)isstring(x)||iscellstr(x)||ischar(x));
p.addParameter('MaxSessions', 3, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('MinAreas', 6, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('MaxAreas', 12, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('MinDuration', 600, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('NumObservations', 2, @(x)isnumeric(x)&&isscalar(x)&&x>=2);
p.addParameter('OverlapFraction', 0.35, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<1);
p.addParameter('NumRepeats', 1, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('Seed', 1, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('Parallel', false, @(x)islogical(x)&&isscalar(x));
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
    opt.OutputFile = fullfile(outDir, 'within_session_subset_stitching.mat');
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
    fprintf('[ibl_validate_within_session_subset_stitching] validating %d sessions\n', numel(sessions));
end

for k = 1:numel(sessions)
    sessionId = string(sessions(k).session_id);
    if cfg.verbose
        fprintf('[subset-stitch validation] %d/%d %s\n', k, numel(sessions), sessionId);
    end

    pcaFile = fullfile(cfg.area_pca_dir, char(sessionId + ".mat"));
    tmp = load(pcaFile, 'summary');
    pcaSummary = tmp.summary;

    areaNames = local_select_areas(pcaSummary, qc.largest_component, opt);
    if numel(areaNames) < opt.MinAreas
        if cfg.verbose
            fprintf('  skip: only %d selected areas\n', numel(areaNames));
        end
        continue;
    end

    regionTable = local_region_table(cfg, qc, pcaSummary);
    areaStruct = ibl_load_area_spikes(pcaSummary.session_root, regionTable, cfg);
    [allNames, ord] = sort(string({areaStruct.area}));
    areaStruct = areaStruct(ord);
    [tf, idxStruct] = ismember(areaNames, allNames);
    assert(all(tf), 'Could not reload all selected areas for %s.', pcaSummary.session_id);
    areaStruct = areaStruct(idxStruct);

    windows = local_observation_windows(pcaSummary.sp_dur, opt.NumObservations);
    [Sfull, f] = local_compute_observation_csd(areaStruct, areaNames, cfg, [0 pcaSummary.sp_dur]);
    SfullMetric = local_normalize_for_metric(Sfull, cfg.stitch.normalize);
    infoflowFullControl = [];
    if opt.ComputePDC
        pdcParams = cfg.stitch;
        pdcParams.verbose = false;
        pdcParams.regularizer = 'eigfloor';
        pdcParams.lambda = 0;
        [~, SregFull] = glasso_precision_estimate(SfullMetric, 0, pdcParams);
        pdcFull = nonparam_pdc_H(SregFull, f, ...
            'maxiter', cfg.chronux.maxiter, 'tol', cfg.chronux.tol);
        infoflowFullControl = ibl_pdc_to_infoflow(pdcFull);
    end

    repeatCells = cell(opt.NumRepeats, 1);
    useParallel = opt.Parallel && opt.NumRepeats > 1;
    if useParallel
        if isempty(gcp('nocreate'))
            parpool;
        end
        parfor r = 1:opt.NumRepeats
            repeatCells{r} = local_run_repeat(r, sessionId, areaNames, areaStruct, ...
                windows, Sfull, SfullMetric, infoflowFullControl, f, cfg, opt);
        end
    else
        for r = 1:opt.NumRepeats
            repeatCells{r} = local_run_repeat(r, sessionId, areaNames, areaStruct, ...
                windows, Sfull, SfullMetric, infoflowFullControl, f, cfg, opt);
        end
    end

    repeatRecords = vertcat(repeatCells{:});
    records = [records; repeatRecords]; %#ok<AGROW>

    if cfg.verbose
        if opt.ComputePDC
            T = struct2table(arrayfun(@(x) x.infoflow_metrics, repeatRecords));
            fprintf('  repeats=%d  infoflow completed corr median=%.3f  IQR=[%.3f %.3f]\n', ...
                opt.NumRepeats, median(T.completed_corr, 'omitnan'), ...
                prctile(T.completed_corr, 25), prctile(T.completed_corr, 75));
        else
            T = struct2table(arrayfun(@(x) x.metrics, repeatRecords));
            fprintf('  repeats=%d  CSD completed corr median=%.3f  IQR=[%.3f %.3f]\n', ...
                opt.NumRepeats, median(T.completed_abs_corr, 'omitnan'), ...
                prctile(T.completed_abs_corr, 25), prctile(T.completed_abs_corr, 75));
        end
        if opt.MakeFigure && ~useParallel
            for r = 1:numel(repeatRecords)
                local_plot_record(repeatRecords(r));
            end
        end
    end
end

val.kind = "within_session_subset_stitching";
val.records = records;
val.summary_table = local_records_to_table(records, opt.ComputePDC);
val.options = opt;
val.created_at = string(datetime('now'));
save(opt.OutputFile, 'val', '-v7.3');
if cfg.verbose
    fprintf('[ibl_validate_within_session_subset_stitching] saved %s\n', opt.OutputFile);
end
end


function T = local_records_to_table(records, hasPdc)
if isempty(records)
    T = table();
    return;
end

n = numel(records);
session_id = strings(n,1);
repeat_id = zeros(n,1);
seed = zeros(n,1);
n_areas = zeros(n,1);
n_observed_pairs = zeros(n,1);
n_completed_pairs = zeros(n,1);
csd_all_abs_corr = nan(n,1);
csd_observed_abs_corr = nan(n,1);
csd_completed_abs_corr = nan(n,1);
csd_all_rel_rmse = nan(n,1);
csd_observed_rel_rmse = nan(n,1);
csd_completed_rel_rmse = nan(n,1);
if hasPdc
    info_all_corr = nan(n,1);
    info_observed_corr = nan(n,1);
    info_completed_corr = nan(n,1);
    info_all_rel_rmse = nan(n,1);
    info_observed_rel_rmse = nan(n,1);
    info_completed_rel_rmse = nan(n,1);
end

for i = 1:n
    session_id(i) = string(records(i).session_id);
    repeat_id(i) = records(i).repeat_id;
    seed(i) = records(i).seed;
    n_areas(i) = numel(records(i).area_names);
    n_observed_pairs(i) = records(i).metrics.n_observed_pairs;
    n_completed_pairs(i) = records(i).metrics.n_completed_pairs;
    csd_all_abs_corr(i) = records(i).metrics.all_abs_corr;
    csd_observed_abs_corr(i) = records(i).metrics.observed_abs_corr;
    csd_completed_abs_corr(i) = records(i).metrics.completed_abs_corr;
    csd_all_rel_rmse(i) = records(i).metrics.all_rel_rmse;
    csd_observed_rel_rmse(i) = records(i).metrics.observed_rel_rmse;
    csd_completed_rel_rmse(i) = records(i).metrics.completed_rel_rmse;
    if hasPdc
        info_all_corr(i) = records(i).infoflow_metrics.all_corr;
        info_observed_corr(i) = records(i).infoflow_metrics.observed_corr;
        info_completed_corr(i) = records(i).infoflow_metrics.completed_corr;
        info_all_rel_rmse(i) = records(i).infoflow_metrics.all_rel_rmse;
        info_observed_rel_rmse(i) = records(i).infoflow_metrics.observed_rel_rmse;
        info_completed_rel_rmse(i) = records(i).infoflow_metrics.completed_rel_rmse;
    end
end

T = table(session_id, repeat_id, seed, n_areas, n_observed_pairs, n_completed_pairs, ...
    csd_all_abs_corr, csd_observed_abs_corr, csd_completed_abs_corr, ...
    csd_all_rel_rmse, csd_observed_rel_rmse, csd_completed_rel_rmse);
if hasPdc
    T.info_all_corr = info_all_corr;
    T.info_observed_corr = info_observed_corr;
    T.info_completed_corr = info_completed_corr;
    T.info_all_rel_rmse = info_all_rel_rmse;
    T.info_observed_rel_rmse = info_observed_rel_rmse;
    T.info_completed_rel_rmse = info_completed_rel_rmse;
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
    pcaFile = fullfile(cfg.area_pca_dir, char(string(sessions(k).session_id) + ".mat"));
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


function regionTable = local_region_table(cfg, qc, pcaSummary)
regionTable = [];
if isfield(cfg, 'region_table')
    regionTable = cfg.region_table;
end
if isempty(regionTable) && isfield(pcaSummary, 'region_table')
    regionTable = pcaSummary.region_table;
end
if isempty(regionTable) && isfield(qc, 'scan') && isfield(qc.scan, 'region_table')
    regionTable = qc.scan.region_table;
end
if isempty(regionTable)
    tmp = load(cfg.scan_file, 'scan');
    regionTable = tmp.scan.region_table;
end
end


function rec = local_run_repeat(repeatId, sessionId, areaNames, areaStruct, ...
    windows, Sfull, SfullMetric, infoflowFullControl, f, cfg, opt)
rng(opt.Seed + repeatId - 1);
subsets = local_area_subsets(numel(areaNames), opt.NumObservations, opt.OverlapFraction);

Sblocks = cell(1, opt.NumObservations);
recset = cell(1, opt.NumObservations);
obsSummaries = struct('window', {}, 'area_names', {}, 'indices', {});
for u = 1:opt.NumObservations
    idx = subsets{u};
    recset{u} = idx(:).';
    Sblocks{u} = local_compute_observation_csd(areaStruct(idx), areaNames(idx), cfg, windows(u,:));
    obsSummaries(u).window = windows(u,:);
    obsSummaries(u).area_names = areaNames(idx);
    obsSummaries(u).indices = idx(:).';
end

stitchParams = cfg.stitch;
stitchParams.verbose = false;
[~, Sstitch, SregStitch, f2, meta] = stitch_spectra_blocks(Sblocks, recset, f, stitchParams);
assert(numel(f) == numel(f2) && max(abs(f(:)-f2(:))) < 1e-10, ...
    'Stitching changed frequency grid for %s.', sessionId);

observedPair = ~meta.missing_mask;
rec = struct();
rec.session_id = sessionId;
rec.repeat_id = repeatId;
rec.seed = opt.Seed + repeatId - 1;
rec.area_names = areaNames;
rec.freqs = f;
rec.observations = obsSummaries;
rec.observed_pair_count = meta.count_mat;
rec.missing_mask = meta.missing_mask;
rec.S_full_control = Sfull;
rec.S_full_control_metric = SfullMetric;
rec.S_stitched = Sstitch;
rec.metric_normalize = string(cfg.stitch.normalize);
rec.metrics = local_csd_metrics(Sstitch, SfullMetric, observedPair);

if opt.ComputePDC
    pdcStitch = nonparam_pdc_H(SregStitch, f, ...
        'maxiter', cfg.chronux.maxiter, 'tol', cfg.chronux.tol);
    rec.infoflow_stitch = ibl_pdc_to_infoflow(pdcStitch);
    rec.infoflow_full_control = infoflowFullControl;
    rec.infoflow_metrics = local_infoflow_metrics( ...
        rec.infoflow_stitch, rec.infoflow_full_control, observedPair);
end
end


function windows = local_observation_windows(spDur, nObs)
edges = linspace(0, spDur, nObs + 1);
windows = [edges(1:end-1).' edges(2:end).'];
end


function subsets = local_area_subsets(nAreas, nObs, overlapFraction)
if nObs == 2
    nOverlap = max(1, min(nAreas - 1, round(overlapFraction * nAreas)));
    ord = randperm(nAreas);
    overlap = ord(1:nOverlap);
    rest = ord(nOverlap+1:end);
    cut = ceil(numel(rest) / 2);
    subsets = cell(1, 2);
    subsets{1} = sort([overlap rest(1:cut)]);
    subsets{2} = sort([overlap rest(cut+1:end)]);
    return;
end

nPerObs = max(2, ceil((1 + overlapFraction) * nAreas / nObs));
subsets = cell(1, nObs);
areaBins = mod(randperm(nAreas) - 1, nObs) + 1;
for u = 1:nObs
    mustHave = find(areaBins == u);
    pool = setdiff(1:nAreas, mustHave);
    addN = max(0, nPerObs - numel(mustHave));
    add = pool(randperm(numel(pool), min(addN, numel(pool))));
    subsets{u} = sort([mustHave add]);
end

covered = unique([subsets{:}]);
missing = setdiff(1:nAreas, covered);
for k = 1:numel(missing)
    u = mod(k - 1, nObs) + 1;
    subsets{u} = sort(unique([subsets{u} missing(k)]));
end
end


function [Sarea, f] = local_compute_observation_csd(areaStruct, areaNames, cfg, win)
Tmax = win(2) - win(1);
nAreas = numel(areaNames);
Sarea = [];
f = [];
loadings = cell(nAreas, 1);

for a = 1:nAreas
    data = local_times_to_struct(areaStruct(a).spike_times, win);
    [Saa, fa] = local_point_spectrum(data, Tmax, cfg.chronux);
    [Saa, fa] = ibl_reduce_frequency_grid(Saa, fa, cfg.target_n_freqs);
    if isempty(f)
        f = fa(:);
        Sarea = zeros(nAreas, nAreas, numel(f));
    else
        assert(numel(f) == numel(fa) && max(abs(f(:)-fa(:))) < 1e-10, ...
            'Frequency mismatch while computing observation auto spectra.');
    end
    [vfix, ~] = ibl_fixed_spectral_loading(Saa, f, cfg.pc_band);
    loadings{a} = vfix;
    Sarea(a,a,:) = ibl_project_cross_spectrum_fixed(Saa, vfix, vfix);
end

for a = 1:nAreas
    for b = (a+1):nAreas
        data = [local_times_to_struct(areaStruct(a).spike_times, win); ...
                local_times_to_struct(areaStruct(b).spike_times, win)];
        [SabFull, fab] = local_point_spectrum(data, Tmax, cfg.chronux);
        [SabFull, fab] = ibl_reduce_frequency_grid(SabFull, fab, cfg.target_n_freqs);
        assert(numel(f) == numel(fab) && max(abs(f(:)-fab(:))) < 1e-10, ...
            'Frequency mismatch while computing observation cross spectra.');
        na = areaStruct(a).n_neurons;
        nb = areaStruct(b).n_neurons;
        Sab = SabFull(1:na, na+1:na+nb, :);
        vab = ibl_project_cross_spectrum_fixed(Sab, loadings{a}, loadings{b});
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


function metrics = local_csd_metrics(Sstitch, Sfull, observedPair)
off = ~eye(size(observedPair));
observed = observedPair & off;
completed = ~observedPair & off;
metrics.all_abs_corr = local_abs_corr(Sstitch, Sfull, off);
metrics.observed_abs_corr = local_abs_corr(Sstitch, Sfull, observed);
metrics.completed_abs_corr = local_abs_corr(Sstitch, Sfull, completed);
metrics.all_rel_rmse = local_rel_rmse(Sstitch, Sfull, off);
metrics.observed_rel_rmse = local_rel_rmse(Sstitch, Sfull, observed);
metrics.completed_rel_rmse = local_rel_rmse(Sstitch, Sfull, completed);
metrics.n_observed_pairs = nnz(triu(observed, 1));
metrics.n_completed_pairs = nnz(triu(completed, 1));
end


function S = local_normalize_for_metric(S, normalizeMode)
mode = lower(string(normalizeMode));
switch mode
    case "none"
        return;
    case {"coherence", "coh"}
        n = size(S, 1);
        for f = 1:size(S, 3)
            A = (S(:,:,f) + S(:,:,f)') / 2;
            p = real(diag(A));
            p(~isfinite(p) | p <= 0) = nan;
            denom = sqrt(p * p.');
            C = A ./ denom;
            C(~isfinite(C)) = 0;
            C(1:n+1:end) = 1;
            S(:,:,f) = (C + C') / 2;
        end
    otherwise
        error('Unknown metric normalize mode: %s', normalizeMode);
end
end


function metrics = local_infoflow_metrics(Fstitch, Ffull, observedPair)
off = ~eye(size(observedPair));
observed = observedPair & off;
completed = ~observedPair & off;
metrics.all_corr = local_vec_corr(Fstitch(off), Ffull(off));
metrics.observed_corr = local_vec_corr(Fstitch(observed), Ffull(observed));
metrics.completed_corr = local_vec_corr(Fstitch(completed), Ffull(completed));
metrics.all_rel_rmse = norm(Fstitch(off) - Ffull(off)) / max(norm(Ffull(off)), eps);
metrics.observed_rel_rmse = norm(Fstitch(observed) - Ffull(observed)) / max(norm(Ffull(observed)), eps);
metrics.completed_rel_rmse = norm(Fstitch(completed) - Ffull(completed)) / max(norm(Ffull(completed)), eps);
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
    va(pos:pos+nf-1) = squeeze(abs(A(i,j,:)));
    vb(pos:pos+nf-1) = squeeze(abs(B(i,j,:)));
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


function local_plot_record(rec)
figure('Color', 'w', 'Position', [100 100 1220 430]);
tiledlayout(1,3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
imagesc(rec.observed_pair_count > 0);
axis image;
title(sprintf('%s subset co-observation', rec.session_id), 'Interpreter', 'none');
set(gca, 'XTick', 1:numel(rec.area_names), 'XTickLabel', rec.area_names, ...
    'YTick', 1:numel(rec.area_names), 'YTickLabel', rec.area_names, ...
    'XTickLabelRotation', 90, 'TickLabelInterpreter', 'none');

nexttile;
bar([rec.metrics.all_abs_corr, rec.metrics.observed_abs_corr, rec.metrics.completed_abs_corr; ...
     rec.metrics.all_rel_rmse, rec.metrics.observed_rel_rmse, rec.metrics.completed_rel_rmse]);
set(gca, 'XTickLabel', {'abs corr', 'rel RMSE'});
legend({'all', 'observed', 'completed'}, 'Location', 'best');
title('Stitched CSD vs full simultaneous control');

nexttile;
if isfield(rec, 'infoflow_metrics')
    bar([rec.infoflow_metrics.all_corr, rec.infoflow_metrics.observed_corr, rec.infoflow_metrics.completed_corr; ...
         rec.infoflow_metrics.all_rel_rmse, rec.infoflow_metrics.observed_rel_rmse, rec.infoflow_metrics.completed_rel_rmse]);
    set(gca, 'XTickLabel', {'corr', 'rel RMSE'});
    legend({'all', 'observed', 'completed'}, 'Location', 'best');
    title('Infoflow vs full simultaneous control');
else
    axis off;
    text(0.1, 0.5, 'PDC/infoflow skipped');
end
end
