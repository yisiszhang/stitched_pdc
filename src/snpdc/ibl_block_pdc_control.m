function ctrl = ibl_block_pdc_control(result, cfg, allenMatrixCsv, varargin)
%IBL_BLOCK_PDC_CONTROL Per-session co-observed PDC control.
%
%   ctrl = ibl_block_pdc_control(result, cfg, allenMatrixCsv)
%
% Computes PDC independently inside each final stitched session block using
% only co-observed areas. Directed information-flow entries are then averaged
% across sessions for pairs observed more than once. Allen comparison is
% restricted to directly observed pairs only.

p = inputParser;
p.addParameter('Band', [1 80], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('RestrictToAllenSupported', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('AllenPositiveThreshold', 0, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('Regularizer', 'eigfloor', @(x)ischar(x)||isstring(x));
p.addParameter('Lambda', 0, @(x)isnumeric(x));
p.addParameter('MakeFigure', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('OutputFile', "", @(x)isstring(x)||ischar(x));
p.addParameter('UseCache', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('ForceRecompute', false, @(x)islogical(x)&&isscalar(x));
p.parse(varargin{:});
opt = p.Results;

areaNames = string(result.area_names(:));
nAreas = numel(areaNames);
freqs = result.freqs(:);
bandMask = freqs >= opt.Band(1) & freqs <= min(opt.Band(2), freqs(end));
assert(any(bandMask), 'No frequencies fall inside [%g %g] Hz.', opt.Band(1), opt.Band(2));

cacheFile = string(opt.OutputFile);
if opt.UseCache && ~opt.ForceRecompute && strlength(cacheFile) > 0 && exist(cacheFile, 'file') == 2
    loaded = load(cacheFile, 'ctrl');
    if isfield(loaded, 'ctrl') && local_cache_matches_result(loaded.ctrl, result)
        ctrl = loaded.ctrl;
        stitchedInfo = ibl_pdc_to_infoflow(result.PDC(:,:,bandMask));
        ctrl.block_control_allen = local_compare_observed_to_allen( ...
            ctrl.infoflow_observed_mean, result, allenMatrixCsv, ctrl.observed_mask, opt, "block_control");
        ctrl.stitched_observed_allen = local_compare_observed_to_allen( ...
            stitchedInfo, result, allenMatrixCsv, ctrl.observed_mask, opt, "stitched_observed");
        ctrl.summary_table = local_summary_table(ctrl.block_control_allen, ctrl.stitched_observed_allen);
        ctrl.options = opt;
        ctrl.cache_status = "loaded_block_pdc_recomputed_allen_only";
        ctrl.updated_at = string(datetime('now'));
        save(cacheFile, 'ctrl', '-v7.3');
        if opt.MakeFigure
            local_plot_control(ctrl);
        end
        if isfield(cfg, 'verbose') && cfg.verbose
            fprintf('[ibl_block_pdc_control] loaded cached block PDC and rescored Allen: %s\n', cacheFile);
        end
        return;
    elseif isfield(cfg, 'verbose') && cfg.verbose
        fprintf('[ibl_block_pdc_control] cache ignored because result areas do not match: %s\n', cacheFile);
    end
end

infoSum = zeros(nAreas);
infoCount = zeros(nAreas);
blockRows = struct([]);

for s = 1:numel(result.kept_sessions)
    sessionId = string(result.kept_sessions(s).session_id);
    specFile = fullfile(cfg.cross_spectra_dir, char(sessionId + ".mat"));
    if exist(specFile, 'file') ~= 2
        warning('ibl_block_pdc_control:MissingSessionSpectrum', ...
            'Missing cross-spectrum file for %s: %s', sessionId, specFile);
        continue;
    end

    tmp = load(specFile, 'summary');
    summary = tmp.summary;
    sessAreas = string(summary.area_names(:));
    [sharedAreas, idxSession, idxGlobal] = intersect(sessAreas, areaNames, 'stable');
    if numel(sharedAreas) < 2
        continue;
    end
    Sblock = summary.cross_spectrum(idxSession, idxSession, :);
    assert(numel(summary.freqs) == numel(freqs) && max(abs(summary.freqs(:) - freqs)) < 1e-10, ...
        'Frequency mismatch for block control session %s.', sessionId);

    params = cfg.stitch;
    params.method = 'none';
    params.regularizer = char(opt.Regularizer);
    params.lambda = opt.Lambda;
    params.verbose = false;
    [~, ~, SregBlock, fBlock] = stitch_spectra_blocks({Sblock}, {1:numel(sharedAreas)}, freqs, params);
    PDCblock = nonparam_pdc_H(SregBlock, fBlock, ...
        'maxiter', cfg.chronux.maxiter, 'tol', cfg.chronux.tol);
    infoBlock = ibl_pdc_to_infoflow(PDCblock(:,:,bandMask));

    off = ~eye(numel(sharedAreas));
    tmpCount = double(off);
    infoSum(idxGlobal, idxGlobal) = infoSum(idxGlobal, idxGlobal) + infoBlock .* tmpCount;
    infoCount(idxGlobal, idxGlobal) = infoCount(idxGlobal, idxGlobal) + tmpCount;

    row.session_id = sessionId;
    row.n_areas = numel(sharedAreas);
    row.areas = strjoin(sharedAreas(:).', "; ");
    row.n_directed_pairs = nnz(off);
    blockRows = local_append_struct(blockRows, row);

    if isfield(cfg, 'verbose') && cfg.verbose
        fprintf('[ibl_block_pdc_control] %3d/%3d %s areas=%d\n', ...
            s, numel(result.kept_sessions), sessionId, numel(sharedAreas));
    end
end

infoMean = infoSum ./ max(infoCount, 1);
infoMean(infoCount == 0) = NaN;
infoMean(1:nAreas+1:end) = NaN;
observedMask = infoCount > 0;
observedMask(1:nAreas+1:end) = false;

allenCmp = local_compare_observed_to_allen(infoMean, result, allenMatrixCsv, observedMask, opt, "block_control");
stitchedInfo = ibl_pdc_to_infoflow(result.PDC(:,:,bandMask));
stitchedCmp = local_compare_observed_to_allen(stitchedInfo, result, allenMatrixCsv, observedMask, opt, "stitched_observed");

ctrl.area_names = areaNames;
ctrl.freqs = freqs;
ctrl.band = opt.Band;
ctrl.infoflow_observed_mean = infoMean;
ctrl.observed_count = infoCount;
ctrl.observed_mask = observedMask;
ctrl.block_table = local_struct_to_table(blockRows);
ctrl.block_control_allen = allenCmp;
ctrl.stitched_observed_allen = stitchedCmp;
ctrl.summary_table = local_summary_table(allenCmp, stitchedCmp);
ctrl.options = opt;
ctrl.cache_status = "computed_block_pdc";
ctrl.created_at = string(datetime('now'));

if strlength(string(opt.OutputFile)) > 0
    outFile = string(opt.OutputFile);
    outDir = fileparts(char(outFile));
    if ~isempty(outDir) && exist(outDir, 'dir') ~= 7
        mkdir(outDir);
    end
    save(outFile, 'ctrl', '-v7.3');
end

if opt.MakeFigure
    local_plot_control(ctrl);
end
end

function tf = local_cache_matches_result(ctrl, result)
tf = isfield(ctrl, 'infoflow_observed_mean') && ...
    isfield(ctrl, 'observed_mask') && ...
    isfield(ctrl, 'area_names') && ...
    isequal(string(ctrl.area_names(:)), string(result.area_names(:))) && ...
    isequal(size(ctrl.infoflow_observed_mean), [numel(result.area_names) numel(result.area_names)]) && ...
    isequal(size(ctrl.observed_mask), [numel(result.area_names) numel(result.area_names)]);
end

function cmp = local_compare_observed_to_allen(infoflow, result, allenMatrixCsv, observedMask, opt, label)
T = readtable(allenMatrixCsv, 'TextType', 'string', 'VariableNamingRule', 'preserve');
allenTargets = string(T{:,1});
allenSources = string(T.Properties.VariableNames(2:end));
allenVals = table2array(T(:,2:end));

resultAreas = string(result.area_names(:));
shared = intersect(resultAreas, intersect(allenTargets, allenSources), 'stable');
[~, idxResult] = ismember(shared, resultAreas);
[~, idxAllenRows] = ismember(shared, allenTargets);
[~, idxAllenCols] = ismember(shared, allenSources);

F = infoflow(idxResult, idxResult);
A = allenVals(idxAllenRows, idxAllenCols);
O = observedMask(idxResult, idxResult);

if opt.RestrictToAllenSupported
    off = ~eye(numel(shared));
    support = isfinite(A) & off;
    keepArea = any(support, 1).' | any(support, 2);
    shared = shared(keepArea);
    F = F(keepArea, keepArea);
    A = A(keepArea, keepArea);
    O = O(keepArea, keepArea);
end

offMask = ~eye(numel(shared));
pairMask = offMask & O & isfinite(F) & isfinite(A);
scores = F(pairMask);
allenWeights = A(pairMask);
metrics = local_binary_edge_metrics(scores, allenWeights, opt.AllenPositiveThreshold);

cmp.label = string(label);
cmp.shared_areas = shared;
cmp.infoflow = F;
cmp.allen = A;
cmp.observed_pair = O;
cmp.valid_pair_mask = pairMask;
cmp.valid_pairs = metrics.n_pairs;
cmp.n_positive = metrics.n_positive;
cmp.n_negative = metrics.n_negative;
cmp.positive_rate = metrics.positive_rate;
cmp.roc_auc = metrics.roc_auc;
cmp.average_precision = metrics.average_precision;
cmp.pr_auc_trapz = metrics.pr_auc_trapz;
cmp.precision_at_npos = metrics.precision_at_npos;
cmp.recall_at_npos = metrics.recall_at_npos;
cmp.roc = struct('fpr', metrics.fpr, 'tpr', metrics.tpr);
cmp.pr = struct('recall', metrics.recall, 'precision', metrics.precision);
cmp.allen_positive_threshold = opt.AllenPositiveThreshold;
cmp.allen_matrix_csv = string(allenMatrixCsv);
end

function T = local_summary_table(blockCmp, stitchedCmp)
method = [blockCmp.label; stitchedCmp.label];
valid_pairs = [blockCmp.valid_pairs; stitchedCmp.valid_pairs];
n_positive = [blockCmp.n_positive; stitchedCmp.n_positive];
n_negative = [blockCmp.n_negative; stitchedCmp.n_negative];
roc_auc = [blockCmp.roc_auc; stitchedCmp.roc_auc];
average_precision = [blockCmp.average_precision; stitchedCmp.average_precision];
precision_at_npos = [blockCmp.precision_at_npos; stitchedCmp.precision_at_npos];
recall_at_npos = [blockCmp.recall_at_npos; stitchedCmp.recall_at_npos];
T = table(method, valid_pairs, n_positive, n_negative, roc_auc, ...
    average_precision, precision_at_npos, recall_at_npos);
end

function metrics = local_binary_edge_metrics(scores, allenWeights, threshold)
scores = scores(:);
allenWeights = allenWeights(:);
labels = allenWeights > threshold;
finiteMask = isfinite(scores) & isfinite(allenWeights);
scores = scores(finiteMask);
labels = labels(finiteMask);

metrics.n_pairs = numel(labels);
metrics.n_positive = nnz(labels);
metrics.n_negative = nnz(~labels);
metrics.positive_rate = metrics.n_positive / max(numel(labels), 1);
metrics.has_both_classes = metrics.n_positive > 0 && metrics.n_negative > 0;
metrics.roc_auc = NaN;
metrics.average_precision = NaN;
metrics.pr_auc_trapz = NaN;
metrics.precision_at_npos = NaN;
metrics.recall_at_npos = NaN;
metrics.fpr = [0; 1];
metrics.tpr = [0; 1];
metrics.recall = [0; 1];
metrics.precision = [metrics.positive_rate; metrics.positive_rate];

if ~metrics.has_both_classes
    return;
end

[~, ord] = sort(scores, 'descend');
labelsSorted = labels(ord);

ranks = local_average_ranks(scores);
rankSumPositive = sum(ranks(labels));
metrics.roc_auc = (rankSumPositive - metrics.n_positive * (metrics.n_positive + 1) / 2) / ...
    (metrics.n_positive * metrics.n_negative);

tp = cumsum(labelsSorted);
fp = cumsum(~labelsSorted);
metrics.tpr = [0; tp / metrics.n_positive];
metrics.fpr = [0; fp / metrics.n_negative];

metrics.recall = [0; tp / metrics.n_positive];
metrics.precision = [1; tp ./ (tp + fp)];
deltaRecall = diff(metrics.recall);
metrics.average_precision = sum(deltaRecall .* metrics.precision(2:end));
metrics.pr_auc_trapz = trapz(metrics.recall, metrics.precision);

topK = min(metrics.n_positive, numel(labelsSorted));
topLabels = labelsSorted(1:topK);
metrics.precision_at_npos = mean(topLabels);
metrics.recall_at_npos = nnz(topLabels) / metrics.n_positive;
end

function ranks = local_average_ranks(x)
[xs, ord] = sort(x(:), 'ascend');
ranksSorted = zeros(size(xs));
i = 1;
while i <= numel(xs)
    j = i;
    while j < numel(xs) && xs(j + 1) == xs(i)
        j = j + 1;
    end
    ranksSorted(i:j) = (i + j) / 2;
    i = j + 1;
end
ranks = zeros(size(xs));
ranks(ord) = ranksSorted;
end

function local_plot_control(ctrl)
figure('Color', 'w', 'Position', [120 120 980 430]);
tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(ctrl.block_control_allen.roc.fpr, ctrl.block_control_allen.roc.tpr, ...
    'LineWidth', 1.8);
hold on;
plot(ctrl.stitched_observed_allen.roc.fpr, ctrl.stitched_observed_allen.roc.tpr, ...
    'LineWidth', 1.8);
plot([0 1], [0 1], 'k--', 'LineWidth', 1);
axis square; grid on;
xlabel('False positive rate');
ylabel('True positive rate');
legend({sprintf('Block PDC AUC %.3f', ctrl.block_control_allen.roc_auc), ...
        sprintf('Stitched observed AUC %.3f', ctrl.stitched_observed_allen.roc_auc), ...
        'chance'}, 'Location', 'southeast');
title('Observed-pair Allen ROC');

nexttile;
imagesc(ctrl.infoflow_observed_mean);
axis image;
colorbar;
set(gca, 'XTick', 1:numel(ctrl.area_names), 'XTickLabel', ctrl.area_names, ...
    'YTick', 1:numel(ctrl.area_names), 'YTickLabel', ctrl.area_names, ...
    'XTickLabelRotation', 90, 'TickLabelInterpreter', 'none', 'FontSize', 5);
xlabel('Source area');
ylabel('Target area');
title('Per-block mean information flow');
end

function rows = local_append_struct(rows, row)
if isempty(rows)
    rows = row;
else
    rows(end+1) = row; %#ok<AGROW>
end
end

function T = local_struct_to_table(rows)
if isempty(rows)
    T = table();
else
    T = struct2table(rows);
end
end
