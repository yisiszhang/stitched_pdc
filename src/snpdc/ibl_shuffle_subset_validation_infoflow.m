function shuffle = ibl_shuffle_subset_validation_infoflow(validationFile, varargin)
%IBL_SHUFFLE_SUBSET_VALIDATION_INFOFLOW Post-hoc shuffle null for validation infoflow.
%
%   shuffle = ibl_shuffle_subset_validation_infoflow(validationFile)
%
% Loads a saved val struct from ibl_validate_within_session_subset_stitching
% and shuffles stitched information-flow values across completed pairs. This
% preserves the full-control matrix and the number of completed edges, but
% destroys area-pair correspondence. No spectra/PDC are recomputed.

p = inputParser;
p.addParameter('NumShuffles', 1000, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('Seed', 20260510, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('OutputFile', "", @(x)isstring(x)||ischar(x));
p.addParameter('MakeFigure', true, @(x)islogical(x)&&isscalar(x));
p.parse(varargin{:});
opt = p.Results;

tmp = load(validationFile, 'val');
val = tmp.val;
assert(isfield(val, 'records') && ~isempty(val.records), ...
    'Validation file does not contain non-empty val.records.');
assert(isfield(val.records(1), 'infoflow_stitch') && isfield(val.records(1), 'infoflow_full_control'), ...
    'Saved validation records do not contain infoflow matrices. Rerun with ComputePDC=true.');

rng(opt.Seed);
nRecords = numel(val.records);
rows = struct('session_id', {}, 'repeat_id', {}, 'seed', {}, ...
    'n_completed_edges', {}, 'observed_corr', {}, 'completed_corr', {}, ...
    'null_mean', {}, 'null_median', {}, 'null_sd', {}, ...
    'null_p95', {}, 'p_right', {}, 'z_score', {});
nullCorr = nan(nRecords, opt.NumShuffles);

for r = 1:nRecords
    rec = val.records(r);
    completedMask = rec.missing_mask & ~eye(size(rec.missing_mask));
    observedMask = ~rec.missing_mask & ~eye(size(rec.missing_mask));
    realCompleted = local_vec_corr(rec.infoflow_stitch(completedMask), ...
        rec.infoflow_full_control(completedMask));
    realObserved = local_vec_corr(rec.infoflow_stitch(observedMask), ...
        rec.infoflow_full_control(observedMask));

    x = rec.infoflow_stitch(completedMask);
    y = rec.infoflow_full_control(completedMask);
    x = real(x(:));
    y = real(y(:));
    good = isfinite(x) & isfinite(y);
    x = x(good);
    y = y(good);
    for s = 1:opt.NumShuffles
        nullCorr(r,s) = local_vec_corr(x(randperm(numel(x))), y);
    end

    row.session_id = string(rec.session_id);
    row.repeat_id = rec.repeat_id;
    row.seed = rec.seed;
    row.n_completed_edges = numel(x);
    row.observed_corr = realObserved;
    row.completed_corr = realCompleted;
    row.null_mean = mean(nullCorr(r,:), 'omitnan');
    row.null_median = median(nullCorr(r,:), 'omitnan');
    row.null_sd = std(nullCorr(r,:), 'omitnan');
    row.null_p95 = prctile(nullCorr(r,:), 95);
    row.p_right = (1 + nnz(nullCorr(r,:) >= realCompleted)) / (1 + opt.NumShuffles);
    row.z_score = (realCompleted - row.null_mean) / max(row.null_sd, eps);
    rows(end+1) = row; %#ok<AGROW>
end

summaryTable = struct2table(rows);
sessionSummary = local_session_summary(summaryTable);

shuffle.validation_file = string(validationFile);
shuffle.options = opt;
shuffle.summary_table = summaryTable;
shuffle.session_summary = sessionSummary;
shuffle.null_completed_corr = nullCorr;
shuffle.created_at = string(datetime('now'));

if strlength(string(opt.OutputFile)) == 0
    [folder, base] = fileparts(validationFile);
    opt.OutputFile = fullfile(folder, char(string(base) + "_shuffle_null.mat"));
end
shuffle.options.OutputFile = opt.OutputFile;
save(opt.OutputFile, 'shuffle', '-v7.3');

if opt.MakeFigure
    local_plot_shuffle(shuffle);
end

fprintf('[ibl_shuffle_subset_validation_infoflow] saved %s\n', opt.OutputFile);
disp(sessionSummary);
end


function T = local_session_summary(summaryTable)
sessions = unique(summaryTable.session_id, 'stable');
rows = struct('session_id', {}, 'n_repeats', {}, ...
    'completed_corr_median', {}, 'completed_corr_iqr_low', {}, 'completed_corr_iqr_high', {}, ...
    'null_median', {}, 'p_right_median', {}, 'frac_p_lt_0p05', {}, 'z_median', {});
for i = 1:numel(sessions)
    mask = summaryTable.session_id == sessions(i);
    vals = summaryTable.completed_corr(mask);
    nullVals = summaryTable.null_median(mask);
    pvals = summaryTable.p_right(mask);
    zvals = summaryTable.z_score(mask);
    row.session_id = sessions(i);
    row.n_repeats = nnz(mask);
    row.completed_corr_median = median(vals, 'omitnan');
    row.completed_corr_iqr_low = prctile(vals, 25);
    row.completed_corr_iqr_high = prctile(vals, 75);
    row.null_median = median(nullVals, 'omitnan');
    row.p_right_median = median(pvals, 'omitnan');
    row.frac_p_lt_0p05 = mean(pvals < 0.05, 'omitnan');
    row.z_median = median(zvals, 'omitnan');
    rows(end+1) = row; %#ok<AGROW>
end
T = struct2table(rows);
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


function local_plot_shuffle(shuffle)
T = shuffle.summary_table;
sessions = unique(T.session_id, 'stable');
figure('Color', 'w', 'Position', [120 120 1050 430]);
tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on;
for i = 1:numel(sessions)
    mask = T.session_id == sessions(i);
    x = i + 0.12 * randn(nnz(mask), 1);
    scatter(x, T.completed_corr(mask), 35, 'filled', 'MarkerFaceAlpha', 0.65);
    plot([i-0.25 i+0.25], repmat(median(T.completed_corr(mask), 'omitnan'), 1, 2), ...
        'k-', 'LineWidth', 2);
    plot([i-0.25 i+0.25], repmat(median(T.null_median(mask), 'omitnan'), 1, 2), ...
        'r--', 'LineWidth', 1.5);
end
set(gca, 'XTick', 1:numel(sessions), 'XTickLabel', sessions, ...
    'XTickLabelRotation', 30, 'TickLabelInterpreter', 'none');
ylabel('Completed-pair infoflow correlation');
title('Real repeats vs shuffle-null median');
grid on;

nexttile;
histogram(T.p_right, 'BinEdges', 0:0.05:1);
xlabel('Right-tail shuffle p-value');
ylabel('Repeat count');
title('Per-repeat completed-pair shuffle test');
grid on;
end
