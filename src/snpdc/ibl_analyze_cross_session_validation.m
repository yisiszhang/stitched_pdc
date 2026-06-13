function out = ibl_analyze_cross_session_validation(val, varargin)
%IBL_ANALYZE_CROSS_SESSION_VALIDATION Relate CSD reproducibility to PDC quality.
%
%   out = ibl_analyze_cross_session_validation(val)
%   out = ibl_analyze_cross_session_validation('/path/to/validation.mat')
%
% Uses saved cross-session validation results. No spectra are recomputed.

p = inputParser;
p.addParameter('MinObservedCsdCorr', 0.35, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('MinCompletedCsdCorr', -inf, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('MinObservedPairFraction', 0.4, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('MaxCompletedPairFraction', 0.7, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('MinAnchorFraction', 0.25, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('MakeFigure', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('OutputFile', "", @(x)isstring(x)||ischar(x));
p.parse(varargin{:});
opt = p.Results;

if ischar(val) || isstring(val)
    tmp = load(val, 'val');
    val = tmp.val;
end
assert(isstruct(val) && isfield(val, 'records'), 'Input must be a validation struct or .mat with val.records.');

T = local_records_to_table(val.records);
rel = local_relationship_table(T);
[T, selected] = local_rank_candidates(T, opt);

out.options = opt;
out.candidate_table = T;
out.relationship_table = rel;
out.selected_table = selected;
out.created_at = string(datetime('now'));

if strlength(string(opt.OutputFile)) > 0
    outFile = string(opt.OutputFile);
    outDir = fileparts(char(outFile));
    if ~isempty(outDir) && exist(outDir, 'dir') ~= 7
        mkdir(outDir);
    end
    save(outFile, 'out', '-v7.3');
    fprintf('[ibl_analyze_cross_session_validation] saved %s\n', char(outFile));
end

fprintf('[ibl_analyze_cross_session_validation] candidates=%d  selected=%d\n', height(T), height(selected));
if ~isempty(rel)
    disp(rel);
end

if opt.MakeFigure
    local_plot(T, opt);
end
end


function T = local_records_to_table(records)
n = numel(records);
target_session_id = strings(n,1);
candidate_index = zeros(n,1);
n_sources = zeros(n,1);
n_areas = zeros(n,1);
n_observed_pairs = zeros(n,1);
n_completed_pairs = zeros(n,1);
n_total_pairs = zeros(n,1);
observed_pair_fraction = nan(n,1);
completed_pair_fraction = nan(n,1);
anchor_fraction = nan(n,1);
source_graph_connected = false(n,1);
min_source_overlap = nan(n,1);
mean_source_jaccard = nan(n,1);
csd_all_abs_corr = nan(n,1);
csd_observed_abs_corr = nan(n,1);
csd_completed_abs_corr = nan(n,1);
csd_all_rel_rmse = nan(n,1);
csd_observed_rel_rmse = nan(n,1);
csd_completed_rel_rmse = nan(n,1);
info_all_corr = nan(n,1);
info_observed_corr = nan(n,1);
info_completed_corr = nan(n,1);
info_all_rel_rmse = nan(n,1);
info_observed_rel_rmse = nan(n,1);
info_completed_rel_rmse = nan(n,1);
source_session_ids = strings(n,1);
covered_areas = strings(n,1);

for i = 1:n
    rec = records(i);
    target_session_id(i) = string(rec.target_session_id);
    if isfield(rec, 'candidate_index')
        candidate_index(i) = rec.candidate_index;
    else
        candidate_index(i) = i;
    end
    n_sources(i) = numel(rec.sources);
    n_areas(i) = numel(rec.covered_area_names);
    n_observed_pairs(i) = rec.metrics.n_observed_pairs;
    n_completed_pairs(i) = rec.metrics.n_completed_pairs;
    n_total_pairs(i) = n_observed_pairs(i) + n_completed_pairs(i);
    observed_pair_fraction(i) = n_observed_pairs(i) / max(n_total_pairs(i), 1);
    completed_pair_fraction(i) = n_completed_pairs(i) / max(n_total_pairs(i), 1);
    csd_all_abs_corr(i) = rec.metrics.all_abs_corr;
    csd_observed_abs_corr(i) = rec.metrics.observed_abs_corr;
    csd_completed_abs_corr(i) = rec.metrics.completed_abs_corr;
    csd_all_rel_rmse(i) = rec.metrics.all_rel_rmse;
    csd_observed_rel_rmse(i) = rec.metrics.observed_rel_rmse;
    csd_completed_rel_rmse(i) = rec.metrics.completed_rel_rmse;
    if isfield(rec, 'infoflow_metrics')
        info_all_corr(i) = rec.infoflow_metrics.all_corr;
        info_observed_corr(i) = rec.infoflow_metrics.observed_corr;
        info_completed_corr(i) = rec.infoflow_metrics.completed_corr;
        info_all_rel_rmse(i) = rec.infoflow_metrics.all_rel_rmse;
        info_observed_rel_rmse(i) = rec.infoflow_metrics.observed_rel_rmse;
        info_completed_rel_rmse(i) = rec.infoflow_metrics.completed_rel_rmse;
    end
    if isfield(rec, 'observed_pair_count')
        countMat = rec.observed_pair_count;
        anchor_fraction(i) = nnz(sum(countMat > 0, 2) >= 2) / max(size(countMat,1), 1);
    end
    [source_graph_connected(i), min_source_overlap(i), mean_source_jaccard(i)] = ...
        local_source_overlap_from_record(rec);
    source_session_ids(i) = local_join_source_ids(rec);
    covered_areas(i) = strjoin(string(rec.covered_area_names(:)).', ", ");
end

T = table(candidate_index, target_session_id, n_sources, n_areas, ...
    n_observed_pairs, n_completed_pairs, n_total_pairs, ...
    observed_pair_fraction, completed_pair_fraction, anchor_fraction, ...
    source_graph_connected, min_source_overlap, mean_source_jaccard, ...
    csd_all_abs_corr, csd_observed_abs_corr, csd_completed_abs_corr, ...
    csd_all_rel_rmse, csd_observed_rel_rmse, csd_completed_rel_rmse, ...
    info_all_corr, info_observed_corr, info_completed_corr, ...
    info_all_rel_rmse, info_observed_rel_rmse, info_completed_rel_rmse, ...
    source_session_ids, covered_areas);
end


function rel = local_relationship_table(T)
x = T.csd_observed_abs_corr;
names = ["info_all_corr"; "info_observed_corr"; "info_completed_corr"; ...
         "csd_completed_abs_corr"];
rho = nan(numel(names),1);
pval = nan(numel(names),1);
n = zeros(numel(names),1);
for i = 1:numel(names)
    y = T.(names(i));
    good = isfinite(x) & isfinite(y);
    n(i) = nnz(good);
    if n(i) >= 3
        [rho(i), pval(i)] = corr(x(good), y(good), 'type', 'Spearman');
    end
end
rel = table(names, n, rho, pval, 'VariableNames', ...
    {'outcome', 'n', 'spearman_rho_vs_csd_observed', 'p_value'});
end


function [T, selected] = local_rank_candidates(T, opt)
T.selection_pass = ...
    T.source_graph_connected & ...
    T.observed_pair_fraction >= opt.MinObservedPairFraction & ...
    T.completed_pair_fraction <= opt.MaxCompletedPairFraction & ...
    T.anchor_fraction >= opt.MinAnchorFraction & ...
    T.csd_observed_abs_corr >= opt.MinObservedCsdCorr & ...
    T.csd_completed_abs_corr >= opt.MinCompletedCsdCorr;

T.reconstruction_score = ...
    2.0 * T.csd_observed_abs_corr + ...
    1.0 * max(T.csd_completed_abs_corr, -0.2) + ...
    0.5 * T.observed_pair_fraction + ...
    0.5 * T.anchor_fraction - ...
    0.5 * T.completed_pair_fraction;

T = sortrows(T, {'selection_pass', 'reconstruction_score'}, {'descend', 'descend'});
selected = T(T.selection_pass, :);
end


function [connected, minOverlap, meanJaccard] = local_source_overlap_from_record(rec)
n = numel(rec.sources);
if n <= 1
    connected = true;
    minOverlap = NaN;
    meanJaccard = NaN;
    return;
end
areas = cell(n,1);
for i = 1:n
    if isfield(rec.sources(i), 'area_names')
        areas{i} = string(rec.sources(i).area_names(:));
    else
        areas{i} = strings(0,1);
    end
end
overlap = zeros(n);
jaccard = nan(n);
for i = 1:n
    for j = (i+1):n
        inter = numel(intersect(areas{i}, areas{j}));
        uni = numel(union(areas{i}, areas{j}));
        overlap(i,j) = inter;
        overlap(j,i) = inter;
        if uni > 0
            jaccard(i,j) = inter / uni;
            jaccard(j,i) = jaccard(i,j);
        end
    end
end
adj = overlap >= 1;
adj(1:n+1:end) = true;
connected = local_is_connected(adj);
vals = overlap(triu(true(n), 1));
jvals = jaccard(triu(true(n), 1));
minOverlap = local_nanmin(vals);
meanJaccard = local_nanmean(jvals);
end


function tf = local_is_connected(adj)
n = size(adj, 1);
seen = false(n,1);
queue = 1;
seen(1) = true;
while ~isempty(queue)
    cur = queue(1);
    queue(1) = [];
    nb = find(adj(cur,:) & ~seen.');
    seen(nb) = true;
    queue = [queue nb]; %#ok<AGROW>
end
tf = all(seen);
end


function ids = local_join_source_ids(rec)
if ~isfield(rec, 'sources') || isempty(rec.sources)
    ids = "";
    return;
end
ids = strings(numel(rec.sources),1);
for i = 1:numel(rec.sources)
    ids(i) = string(rec.sources(i).session_id);
end
ids = strjoin(ids.', ", ");
end


function y = local_nanmin(x)
x = x(isfinite(x));
if isempty(x)
    y = NaN;
else
    y = min(x);
end
end


function y = local_nanmean(x)
x = x(isfinite(x));
if isempty(x)
    y = NaN;
else
    y = mean(x);
end
end


function local_plot(T, opt)
figure('Color', 'w', 'Position', [100 100 1220 430]);
tiledlayout(1,3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
scatter(T.csd_observed_abs_corr, T.info_all_corr, 55, T.observed_pair_fraction, 'filled');
xline(opt.MinObservedCsdCorr, 'r-', 'LineWidth', 1.2);
xlabel('Observed-pair CSD abs corr');
ylabel('Infoflow all corr');
title('CSD reproducibility vs all infoflow');
colorbar;
grid on;

nexttile;
scatter(T.csd_observed_abs_corr, T.info_observed_corr, 55, T.anchor_fraction, 'filled');
xline(opt.MinObservedCsdCorr, 'r-', 'LineWidth', 1.2);
xlabel('Observed-pair CSD abs corr');
ylabel('Infoflow observed corr');
title('CSD reproducibility vs observed infoflow');
colorbar;
grid on;

nexttile;
scatter(T.observed_pair_fraction, T.csd_completed_abs_corr, 55, double(T.selection_pass), 'filled');
xline(opt.MinObservedPairFraction, 'r-', 'LineWidth', 1.2);
xlabel('Observed pair fraction');
ylabel('Completed-pair CSD abs corr');
title('Coverage geometry vs completion');
colorbar;
grid on;
end
