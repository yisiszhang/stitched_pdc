function summaryTable = ibl_compare_allen_tracing_sweep(result, allenDir, varargin)
%IBL_COMPARE_ALLEN_TRACING_SWEEP Compare stitched info flow to multiple Allen matrices.
%
%   summaryTable = ibl_compare_allen_tracing_sweep(result, allenDir)

files = dir(fullfile(allenDir, 'allen_tracing_matrix__*.csv'));
assert(~isempty(files), 'No Allen matrix files found in %s', allenDir);

rows = struct('file', {}, 'metric', {}, 'aggregate', {}, ...
    'n_shared_areas', {}, 'n_valid_pairs', {}, ...
    'n_positive', {}, 'n_negative', {}, ...
    'pearson_r', {}, 'spearman_r', {}, ...
    'roc_auc', {}, 'average_precision', {}, ...
    'precision_at_npos', {}, 'recall_at_npos', {}, ...
    'observed_valid_pairs', {}, 'observed_n_positive', {}, ...
    'observed_roc_auc', {}, 'observed_average_precision', {}, ...
    'observed_precision_at_npos', {}, 'observed_recall_at_npos', {}, ...
    'never_observed_valid_pairs', {}, 'never_observed_n_positive', {}, ...
    'never_observed_roc_auc', {}, 'never_observed_average_precision', {}, ...
    'never_observed_precision_at_npos', {}, 'never_observed_recall_at_npos', {});

for k = 1:numel(files)
    file = fullfile(files(k).folder, files(k).name);
    cmp = ibl_compare_allen_tracing(result, file, varargin{:});

    tokens = regexp(files(k).name, '^allen_tracing_matrix__(.+)__(.+)\.csv$', 'tokens', 'once');
    if isempty(tokens)
        metric = "";
        agg = "";
    else
        metric = string(tokens{1});
        agg = string(tokens{2});
    end

    rows(end+1).file = string(files(k).name); %#ok<AGROW>
    rows(end).metric = metric;
    rows(end).aggregate = agg;
    rows(end).n_shared_areas = numel(cmp.shared_areas);
    rows(end).n_valid_pairs = cmp.valid_pairs;
    rows(end).n_positive = cmp.n_positive;
    rows(end).n_negative = cmp.n_negative;
    rows(end).pearson_r = cmp.pearson_r;
    rows(end).spearman_r = cmp.spearman_r;
    rows(end).roc_auc = cmp.roc_auc;
    rows(end).average_precision = cmp.average_precision;
    rows(end).precision_at_npos = cmp.precision_at_npos;
    rows(end).recall_at_npos = cmp.recall_at_npos;
    rows(end).observed_valid_pairs = cmp.observed_valid_pairs;
    rows(end).observed_n_positive = cmp.observed_n_positive;
    rows(end).observed_roc_auc = cmp.observed_roc_auc;
    rows(end).observed_average_precision = cmp.observed_average_precision;
    rows(end).observed_precision_at_npos = cmp.observed_precision_at_npos;
    rows(end).observed_recall_at_npos = cmp.observed_recall_at_npos;
    rows(end).never_observed_valid_pairs = cmp.never_observed_valid_pairs;
    rows(end).never_observed_n_positive = cmp.never_observed_n_positive;
    rows(end).never_observed_roc_auc = cmp.never_observed_roc_auc;
    rows(end).never_observed_average_precision = cmp.never_observed_average_precision;
    rows(end).never_observed_precision_at_npos = cmp.never_observed_precision_at_npos;
    rows(end).never_observed_recall_at_npos = cmp.never_observed_recall_at_npos;
end

summaryTable = struct2table(rows);
summaryTable = sortrows(summaryTable, ...
    {'roc_auc', 'average_precision', 'spearman_r', 'pearson_r'}, ...
    {'descend', 'descend', 'descend', 'descend'});
end
