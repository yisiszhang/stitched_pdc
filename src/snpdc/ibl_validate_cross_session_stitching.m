function val = ibl_validate_cross_session_stitching(crossCand, cfg, varargin)
%IBL_VALIDATE_CROSS_SESSION_STITCHING Validate stitching across animals/sessions.
%
%   val = ibl_validate_cross_session_stitching(crossCand, cfg)
%
% Each candidate has one held-out target session used only as a full
% simultaneous control. Source sessions from other animals are stitched and
% compared against the target on the covered target-area subset.

assert(exist('CrossSpecMatpt', 'file') == 2, ...
    'Chronux CrossSpecMatpt not found. Add Chronux to the MATLAB path.');

p = inputParser;
p.addParameter('CandidateIndex', 1, @(x)isnumeric(x)&&all(x>=1));
p.addParameter('ComputePDC', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('MakeFigure', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('OutputFile', "", @(x)isstring(x)||ischar(x));
p.parse(varargin{:});
opt = p.Results;

if isfield(crossCand, 'sets')
    candidateSets = crossCand.sets;
else
    candidateSets = crossCand;
end
idxList = opt.CandidateIndex(:).';
idxList = idxList(idxList >= 1 & idxList <= numel(candidateSets));
assert(~isempty(idxList), 'No valid CandidateIndex values.');

if strlength(string(opt.OutputFile)) == 0
    outDir = fullfile(cfg.output_root, 'validation');
    if exist(outDir, 'dir') ~= 7
        mkdir(outDir);
    end
    opt.OutputFile = fullfile(outDir, 'cross_session_stitching_validation.mat');
end
outDir = fileparts(char(opt.OutputFile));
if ~isempty(outDir) && exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end

records = struct([]);
for c = idxList
    cand = candidateSets(c);
    if cfg.verbose
        fprintf('[cross-session validation] candidate %d  target=%s  sources=%d\n', ...
            c, cand.target.session_id, numel(cand.sources));
    end
    rec = local_validate_one_candidate(c, cand, cfg, opt);
    records = [records; rec]; %#ok<AGROW>
    if cfg.verbose
        fprintf('  CSD abs corr all/observed/completed = %.3f / %.3f / %.3f\n', ...
            rec.metrics.all_abs_corr, rec.metrics.observed_abs_corr, rec.metrics.completed_abs_corr);
        fprintf('  CSD relRMSE all/observed/completed = %.3f / %.3f / %.3f\n', ...
            rec.metrics.all_rel_rmse, rec.metrics.observed_rel_rmse, rec.metrics.completed_rel_rmse);
        if isfield(rec, 'infoflow_metrics')
            fprintf('  infoflow corr all/observed/completed = %.3f / %.3f / %.3f\n', ...
                rec.infoflow_metrics.all_corr, rec.infoflow_metrics.observed_corr, ...
                rec.infoflow_metrics.completed_corr);
            fprintf('  infoflow relRMSE all/observed/completed = %.3f / %.3f / %.3f\n', ...
                rec.infoflow_metrics.all_rel_rmse, rec.infoflow_metrics.observed_rel_rmse, ...
                rec.infoflow_metrics.completed_rel_rmse);
        end
    end
    if opt.MakeFigure
        local_plot_record(rec);
    end
end

val.kind = "cross_session_stitching";
val.records = records;
val.summary_table = local_records_to_table(records, opt.ComputePDC);
val.options = opt;
val.created_at = string(datetime('now'));
save(opt.OutputFile, 'val', '-v7.3');
if cfg.verbose
    fprintf('[ibl_validate_cross_session_stitching] saved %s\n', opt.OutputFile);
end
end


function rec = local_validate_one_candidate(candidateIndex, cand, cfg, opt)
targetPca = local_load_pca(cand.target.pca_file);
regionTable = local_region_table(cfg, targetPca);

coveredAreas = cand.covered_area_names(:);
targetAreas = cand.target_area_names(:);
[~, coveredIdxInTarget] = ismember(coveredAreas, targetAreas);
coveredIdxInTarget = coveredIdxInTarget(:).';

targetAreaStruct = local_load_selected_area_struct(targetPca, regionTable, cfg, coveredAreas);
[Starget, f] = local_compute_session_area_csd(targetAreaStruct, coveredAreas, cfg, targetPca.sp_dur);

Sblocks = cell(1, numel(cand.sources));
recset = cell(1, numel(cand.sources));
sourceSummariesAll = struct('session_id', {}, 'lab', {}, 'subject', {}, ...
    'area_names', {}, 'indices', {});
for s = 1:numel(cand.sources)
    src = cand.sources(s);
    srcPca = local_load_pca(src.pca_file);
    srcAreas = intersect(string(src.area_names(:)), coveredAreas, 'stable');
    [~, idxInCovered] = ismember(srcAreas, coveredAreas);
    idxInCovered = idxInCovered(:).';
    if numel(srcAreas) < 2
        continue;
    end
    srcAreaStruct = local_load_selected_area_struct(srcPca, regionTable, cfg, srcAreas);
    Sblocks{s} = local_compute_session_area_csd(srcAreaStruct, srcAreas, cfg, srcPca.sp_dur);
    recset{s} = idxInCovered;
    srcRec.session_id = src.session_id;
    srcRec.lab = src.lab;
    srcRec.subject = src.subject;
    srcRec.area_names = srcAreas;
    srcRec.indices = idxInCovered;
    sourceSummariesAll(end+1) = srcRec; %#ok<AGROW>
    if cfg.verbose
        fprintf('  source %d/%d  %s  areas=%d\n', ...
            s, numel(cand.sources), src.session_id, numel(srcAreas));
    end
end

keepBlock = ~cellfun(@isempty, Sblocks);
Sblocks = Sblocks(keepBlock);
recset = recset(keepBlock);
sourceSummaries = sourceSummariesAll;
assert(numel(Sblocks) >= 2, 'Candidate %d has fewer than two usable source blocks.', candidateIndex);

stitchParams = cfg.stitch;
stitchParams.verbose = false;
[~, Sstitch, SregStitch, f2, meta] = stitch_spectra_blocks(Sblocks, recset, f, stitchParams);
assert(numel(f) == numel(f2) && max(abs(f(:)-f2(:))) < 1e-10, ...
    'Stitching changed frequency grid for candidate %d.', candidateIndex);

StargetMetric = local_normalize_for_metric(Starget, cfg.stitch.normalize);
observedPair = ~meta.missing_mask;

rec = struct();
rec.candidate_index = candidateIndex;
rec.target_session_id = cand.target.session_id;
rec.target_lab = cand.target.lab;
rec.target_subject = cand.target.subject;
rec.target_area_names = targetAreas;
rec.covered_area_names = coveredAreas;
rec.covered_idx_in_target = coveredIdxInTarget;
rec.sources = sourceSummaries;
rec.freqs = f;
rec.observed_pair_count = meta.count_mat;
rec.missing_mask = meta.missing_mask;
rec.S_target_full_control = Starget;
rec.S_target_full_control_metric = StargetMetric;
rec.S_stitched_sources = Sstitch;
rec.metric_normalize = string(cfg.stitch.normalize);
rec.metrics = local_csd_metrics(Sstitch, StargetMetric, observedPair);

if opt.ComputePDC
    pdcParams = stitchParams;
    pdcParams.regularizer = 'eigfloor';
    pdcParams.lambda = 0;
    [~, SregTarget] = glasso_precision_estimate(StargetMetric, 0, pdcParams);
    pdcTarget = nonparam_pdc_H(SregTarget, f, ...
        'maxiter', cfg.chronux.maxiter, 'tol', cfg.chronux.tol);
    pdcStitch = nonparam_pdc_H(SregStitch, f, ...
        'maxiter', cfg.chronux.maxiter, 'tol', cfg.chronux.tol);
    rec.infoflow_target_full_control = ibl_pdc_to_infoflow(pdcTarget);
    rec.infoflow_stitched_sources = ibl_pdc_to_infoflow(pdcStitch);
    rec.infoflow_metrics = local_infoflow_metrics( ...
        rec.infoflow_stitched_sources, rec.infoflow_target_full_control, observedPair);
end
end


function summary = local_load_pca(pcaFile)
tmp = load(pcaFile, 'summary');
summary = tmp.summary;
end


function regionTable = local_region_table(cfg, pcaSummary)
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
end


function areaStruct = local_load_selected_area_struct(pcaSummary, regionTable, cfg, areaNames)
allStruct = ibl_load_area_spikes(pcaSummary.session_root, regionTable, cfg);
[allNames, ord] = sort(string({allStruct.area}));
allStruct = allStruct(ord);
[tf, idx] = ismember(areaNames, allNames);
assert(all(tf), 'Could not reload all selected areas for %s.', pcaSummary.session_id);
areaStruct = allStruct(idx);
end


function [Sarea, f] = local_compute_session_area_csd(areaStruct, areaNames, cfg, spDur)
Tmax = spDur;
nAreas = numel(areaNames);
Sarea = [];
f = [];
loadings = cell(nAreas, 1);

for a = 1:nAreas
    data = local_times_to_struct(areaStruct(a).spike_times);
    [Saa, fa] = local_point_spectrum(data, Tmax, cfg.chronux);
    [Saa, fa] = ibl_reduce_frequency_grid(Saa, fa, cfg.target_n_freqs);
    if isempty(f)
        f = fa(:);
        Sarea = zeros(nAreas, nAreas, numel(f));
    else
        assert(numel(f) == numel(fa) && max(abs(f(:)-fa(:))) < 1e-10, ...
            'Frequency mismatch while computing area auto spectra.');
    end
    [vfix, ~] = ibl_fixed_spectral_loading(Saa, f, cfg.pc_band);
    loadings{a} = vfix;
    Sarea(a,a,:) = ibl_project_cross_spectrum_fixed(Saa, vfix, vfix);
end

for a = 1:nAreas
    for b = (a+1):nAreas
        data = [local_times_to_struct(areaStruct(a).spike_times); ...
                local_times_to_struct(areaStruct(b).spike_times)];
        [SabFull, fab] = local_point_spectrum(data, Tmax, cfg.chronux);
        [SabFull, fab] = ibl_reduce_frequency_grid(SabFull, fab, cfg.target_n_freqs);
        assert(numel(f) == numel(fab) && max(abs(f(:)-fab(:))) < 1e-10, ...
            'Frequency mismatch while computing cross-area spectra.');
        na = areaStruct(a).n_neurons;
        nb = areaStruct(b).n_neurons;
        Sab = SabFull(1:na, na+1:na+nb, :);
        vab = ibl_project_cross_spectrum_fixed(Sab, loadings{a}, loadings{b});
        Sarea(a,b,:) = vab;
        Sarea(b,a,:) = conj(vab);
    end
end
end


function data = local_times_to_struct(timesCell)
n = numel(timesCell);
data = struct('times', cell(n,1));
for i = 1:n
    data(i).times = timesCell{i}(:);
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


function T = local_records_to_table(records, hasPdc)
if isempty(records)
    T = table();
    return;
end
n = numel(records);
target_session_id = strings(n,1);
n_sources = zeros(n,1);
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
    target_session_id(i) = records(i).target_session_id;
    n_sources(i) = numel(records(i).sources);
    n_areas(i) = numel(records(i).covered_area_names);
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
T = table(target_session_id, n_sources, n_areas, n_observed_pairs, n_completed_pairs, ...
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


function local_plot_record(rec)
figure('Color', 'w', 'Position', [100 100 1220 430]);
tiledlayout(1,3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
imagesc(rec.observed_pair_count > 0);
axis image;
title(sprintf('Cross-session support: %s', rec.target_session_id), 'Interpreter', 'none');
set(gca, 'XTick', 1:numel(rec.covered_area_names), 'XTickLabel', rec.covered_area_names, ...
    'YTick', 1:numel(rec.covered_area_names), 'YTickLabel', rec.covered_area_names, ...
    'XTickLabelRotation', 90, 'TickLabelInterpreter', 'none');

nexttile;
bar([rec.metrics.all_abs_corr, rec.metrics.observed_abs_corr, rec.metrics.completed_abs_corr; ...
     rec.metrics.all_rel_rmse, rec.metrics.observed_rel_rmse, rec.metrics.completed_rel_rmse]);
set(gca, 'XTickLabel', {'abs corr', 'rel RMSE'});
legend({'all', 'observed', 'completed'}, 'Location', 'best');
title('Source-stitched CSD vs target full control');

nexttile;
if isfield(rec, 'infoflow_metrics')
    bar([rec.infoflow_metrics.all_corr, rec.infoflow_metrics.observed_corr, rec.infoflow_metrics.completed_corr; ...
         rec.infoflow_metrics.all_rel_rmse, rec.infoflow_metrics.observed_rel_rmse, rec.infoflow_metrics.completed_rel_rmse]);
    set(gca, 'XTickLabel', {'corr', 'rel RMSE'});
    legend({'all', 'observed', 'completed'}, 'Location', 'best');
    title('Source-stitched infoflow vs target full control');
else
    axis off;
    text(0.1, 0.5, 'PDC/infoflow skipped');
end
end
