function cmp = ibl_compare_allen_tracing(result, allenMatrixCsv, varargin)
%IBL_COMPARE_ALLEN_TRACING Compare stitched information flow with Allen tracing.
%
%   cmp = ibl_compare_allen_tracing(result, allenMatrixCsv)
%   cmp = ibl_compare_allen_tracing(result, allenMatrixCsv, 'Band', [1 80], 'Order', 'cluster')
%
% Inputs
%   result          stitched result struct from ibl_stitch_saved_spectra
%   allenMatrixCsv  CSV from scripts/fetch_allen_tracing.py
%
% Name-value
%   'Band'   frequency band for PDC integration, default [1 80]
%   'Order'  'alphabetical' or 'cluster', default 'cluster'
%   'RestrictToAllenSupported' logical, default true. Keep only areas with
%                              at least one finite incoming or outgoing
%                              Allen tracing value in the shared matrix.
%   'AllenPositiveThreshold' threshold for binarizing Allen edges, default 0
%   'ObservationMask' optional logical full result-space mask; true means
%                     area pair was directly co-observed in at least one
%                     session. Default inferred from result.
%   'MakeMatrixFigure' logical, default true
%   'MakeRocFigure' logical, default true

p = inputParser;
p.addParameter('Band', [1 80], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('Order', 'cluster', @(x)ischar(x)||isstring(x));
p.addParameter('RestrictToAllenSupported', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('AllenPositiveThreshold', 0, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('ObservationMask', [], @(x)islogical(x)||isempty(x));
p.addParameter('MakeMatrixFigure', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('MakeRocFigure', true, @(x)islogical(x)&&isscalar(x));
p.parse(varargin{:});
band = p.Results.Band;
orderMode = lower(string(p.Results.Order));
restrictToAllen = p.Results.RestrictToAllenSupported;
allenPositiveThreshold = p.Results.AllenPositiveThreshold;
observationMaskFull = p.Results.ObservationMask;
makeMatrixFigure = p.Results.MakeMatrixFigure;
makeRocFigure = p.Results.MakeRocFigure;

T = readtable(allenMatrixCsv, 'TextType', 'string', 'VariableNamingRule', 'preserve');
allenTargets = string(T{:,1});
allenSources = string(T.Properties.VariableNames(2:end));
allenVals = table2array(T(:,2:end));

shared = intersect(string(result.area_names(:)), intersect(allenTargets, allenSources), 'stable');
assert(numel(shared) >= 2, 'Too few shared areas between stitched result and Allen matrix.');

[~, idxResult] = ismember(shared, string(result.area_names(:)));
[~, idxAllenRows] = ismember(shared, allenTargets);
[~, idxAllenCols] = ismember(shared, allenSources);

fMask = result.freqs >= band(1) & result.freqs <= min(band(2), result.freqs(end));
assert(any(fMask), 'No stitched frequencies fall inside [%g %g] Hz.', band(1), band(2));
infoflow = ibl_pdc_to_infoflow(result.PDC(:,:,fMask));
infoflow = infoflow(idxResult, idxResult);
allen = allenVals(idxAllenRows, idxAllenCols);
if isempty(observationMaskFull)
    observationMaskFull = local_observation_mask(result);
end
assert(isequal(size(observationMaskFull), [numel(result.area_names) numel(result.area_names)]), ...
    'ObservationMask must be a logical matrix matching result.area_names.');
observedPair = observationMaskFull(idxResult, idxResult);

if restrictToAllen
    offdiagMask = ~eye(numel(shared));
    allenSupport = isfinite(allen) & offdiagMask;
    keepArea = any(allenSupport, 1).' | any(allenSupport, 2);
    shared = shared(keepArea);
    infoflow = infoflow(keepArea, keepArea);
    allen = allen(keepArea, keepArea);
    observedPair = observedPair(keepArea, keepArea);
end

assert(numel(shared) >= 2, 'Too few Allen-supported shared areas remain for comparison.');

perm = 1:numel(shared);
resolvedOrder = char(orderMode);
switch orderMode
    case "cluster"
        X = infoflow;
        X(1:size(X,1)+1:end) = 0;
        rowNorm = sqrt(sum(X.^2, 2));
        rowNorm(rowNorm == 0) = 1;
        X = X ./ rowNorm;
        sim = X * X.';
        sim = max(min(sim, 1), -1);
        D = 1 - sim;
        D = (D + D.') / 2;
        D(1:size(D,1)+1:end) = 0;
        perm = symamd(sparse(D + D.'));
        resolvedOrder = 'symamd';
    otherwise
        [~, perm] = sort(shared);
        resolvedOrder = 'alphabetical';
end

sharedPlot = shared(perm);
infoflowPlot = infoflow(perm, perm);
allenPlot = allen(perm, perm);
observedPairPlot = observedPair(perm, perm);

offMask = ~eye(numel(sharedPlot));
stitchedVec = infoflowPlot(offMask);
allenVec = allenPlot(offMask);
observedVec = observedPairPlot(offMask);
valid = isfinite(stitchedVec) & isfinite(allenVec);

if nnz(valid) >= 3
    pearsonR = corr(stitchedVec(valid), allenVec(valid), 'type', 'Pearson');
    spearmanR = corr(stitchedVec(valid), allenVec(valid), 'type', 'Spearman');
else
    pearsonR = NaN;
    spearmanR = NaN;
end

metrics = local_binary_edge_metrics(stitchedVec(valid), allenVec(valid), allenPositiveThreshold);
metricsObserved = local_binary_edge_metrics( ...
    stitchedVec(valid & observedVec), allenVec(valid & observedVec), allenPositiveThreshold);
metricsNeverObserved = local_binary_edge_metrics( ...
    stitchedVec(valid & ~observedVec), allenVec(valid & ~observedVec), allenPositiveThreshold);

if makeMatrixFigure
    figure('Color', 'w', 'Position', [90 90 1180 520]);
    tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    imagesc(infoflowPlot);
    axis image;
    colorbar;
    set(gca, 'XTick', 1:numel(sharedPlot), 'XTickLabel', sharedPlot, ...
        'YTick', 1:numel(sharedPlot), 'YTickLabel', sharedPlot, ...
        'XTickLabelRotation', 90, 'TickLabelInterpreter', 'none');
    xlabel('Source area');
    ylabel('Target area');
    title(sprintf('Stitched info flow (%g-%g Hz)', band(1), band(2)));

    nexttile;
    imagesc(allenPlot);
    axis image;
    colorbar;
    set(gca, 'XTick', 1:numel(sharedPlot), 'XTickLabel', sharedPlot, ...
        'YTick', 1:numel(sharedPlot), 'YTickLabel', sharedPlot, ...
        'XTickLabelRotation', 90, 'TickLabelInterpreter', 'none');
    xlabel('Source area');
    ylabel('Target area');
    title('Allen tracing (normalized projection volume)');

    sgtitle(sprintf(['Shared areas=%d, order=%s, Pearson=%.3f, Spearman=%.3f, ', ...
        'ROC-AUC all/obs/never=%.3f/%.3f/%.3f'], ...
        numel(sharedPlot), resolvedOrder, pearsonR, spearmanR, ...
        metrics.roc_auc, metricsObserved.roc_auc, metricsNeverObserved.roc_auc));
end

if makeRocFigure && metrics.has_both_classes
    figure('Color', 'w', 'Position', [120 120 980 420]);
    tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    plot(metrics.fpr, metrics.tpr, 'LineWidth', 1.8);
    hold on;
    if metricsObserved.has_both_classes
        plot(metricsObserved.fpr, metricsObserved.tpr, 'LineWidth', 1.5);
    end
    if metricsNeverObserved.has_both_classes
        plot(metricsNeverObserved.fpr, metricsNeverObserved.tpr, 'LineWidth', 1.5);
    end
    plot([0 1], [0 1], 'k--', 'LineWidth', 1);
    axis square;
    grid on;
    xlabel('False positive rate');
    ylabel('True positive rate');
    title(sprintf('ROC-AUC = %.3f', metrics.roc_auc));
    legend(local_curve_labels('ROC-AUC', metrics, metricsObserved, metricsNeverObserved), ...
        'Location', 'southeast', 'Interpreter', 'none');

    nexttile;
    plot(metrics.recall, metrics.precision, 'LineWidth', 1.8);
    hold on;
    if metricsObserved.has_both_classes
        plot(metricsObserved.recall, metricsObserved.precision, 'LineWidth', 1.5);
    end
    if metricsNeverObserved.has_both_classes
        plot(metricsNeverObserved.recall, metricsNeverObserved.precision, 'LineWidth', 1.5);
    end
    yline(metrics.positive_rate, 'k--', 'LineWidth', 1);
    axis square;
    grid on;
    xlabel('Recall');
    ylabel('Precision');
    title(sprintf('Average precision = %.3f', metrics.average_precision));
    legend(local_curve_labels('AP', metrics, metricsObserved, metricsNeverObserved), ...
        'Location', 'northeast', 'Interpreter', 'none');

    sgtitle(sprintf('Allen-positive edges: Allen > %.4g, n+ = %d, n- = %d', ...
        allenPositiveThreshold, metrics.n_positive, metrics.n_negative));
end

cmp.shared_areas = shared;
cmp.shared_areas_plot = sharedPlot;
cmp.infoflow = infoflow;
cmp.allen = allen;
cmp.infoflow_plot = infoflowPlot;
cmp.allen_plot = allenPlot;
cmp.observed_pair = observedPair;
cmp.observed_pair_plot = observedPairPlot;
cmp.permutation = perm(:);
cmp.band = band;
cmp.pearson_r = pearsonR;
cmp.spearman_r = spearmanR;
cmp.valid_pairs = nnz(valid);
cmp.allen_positive_threshold = allenPositiveThreshold;
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
cmp.observed = local_export_metrics(metricsObserved);
cmp.never_observed = local_export_metrics(metricsNeverObserved);
cmp.observed_valid_pairs = metricsObserved.n_pairs;
cmp.never_observed_valid_pairs = metricsNeverObserved.n_pairs;
cmp.observed_n_positive = metricsObserved.n_positive;
cmp.observed_n_negative = metricsObserved.n_negative;
cmp.observed_roc_auc = metricsObserved.roc_auc;
cmp.observed_average_precision = metricsObserved.average_precision;
cmp.observed_precision_at_npos = metricsObserved.precision_at_npos;
cmp.observed_recall_at_npos = metricsObserved.recall_at_npos;
cmp.never_observed_n_positive = metricsNeverObserved.n_positive;
cmp.never_observed_n_negative = metricsNeverObserved.n_negative;
cmp.never_observed_roc_auc = metricsNeverObserved.roc_auc;
cmp.never_observed_average_precision = metricsNeverObserved.average_precision;
cmp.never_observed_precision_at_npos = metricsNeverObserved.precision_at_npos;
cmp.never_observed_recall_at_npos = metricsNeverObserved.recall_at_npos;
cmp.order_mode = resolvedOrder;
end

function observedPair = local_observation_mask(result)
nAreas = numel(result.area_names);
if isfield(result, 'meta') && isfield(result.meta, 'missing_mask') && ...
        isequal(size(result.meta.missing_mask), [nAreas nAreas])
    observedPair = ~logical(result.meta.missing_mask);
elseif isfield(result, 'recset')
    observedPair = meacount_mat(result.recset, nAreas) > 0;
else
    warning('ibl_compare_allen_tracing:MissingObservationMask', ...
        'Could not infer direct co-observation mask; treating all pairs as observed.');
    observedPair = true(nAreas, nAreas);
end
observedPair(1:nAreas+1:end) = false;
end

function out = local_export_metrics(metrics)
out = struct();
out.n_pairs = metrics.n_pairs;
out.n_positive = metrics.n_positive;
out.n_negative = metrics.n_negative;
out.positive_rate = metrics.positive_rate;
out.roc_auc = metrics.roc_auc;
out.average_precision = metrics.average_precision;
out.pr_auc_trapz = metrics.pr_auc_trapz;
out.precision_at_npos = metrics.precision_at_npos;
out.recall_at_npos = metrics.recall_at_npos;
out.roc = struct('fpr', metrics.fpr, 'tpr', metrics.tpr);
out.pr = struct('recall', metrics.recall, 'precision', metrics.precision);
end

function labels = local_curve_labels(metricName, metrics, metricsObserved, metricsNeverObserved)
if strcmp(metricName, 'AP')
    labels = {sprintf('all %.3f', metrics.average_precision)};
    if metricsObserved.has_both_classes
        labels{end+1} = sprintf('observed %.3f', metricsObserved.average_precision); %#ok<AGROW>
    end
    if metricsNeverObserved.has_both_classes
        labels{end+1} = sprintf('never observed %.3f', metricsNeverObserved.average_precision); %#ok<AGROW>
    end
else
    labels = {sprintf('all %.3f', metrics.roc_auc)};
    if metricsObserved.has_both_classes
        labels{end+1} = sprintf('observed %.3f', metricsObserved.roc_auc); %#ok<AGROW>
    end
    if metricsNeverObserved.has_both_classes
        labels{end+1} = sprintf('never observed %.3f', metricsNeverObserved.roc_auc); %#ok<AGROW>
    end
end
labels{end+1} = 'chance';
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
