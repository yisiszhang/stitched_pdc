function thresholdFit = ibl_fit_csd_threshold_from_validation(T, outcomeName, fallbackThreshold)
%IBL_FIT_CSD_THRESHOLD_FROM_VALIDATION Fit CSD reproducibility threshold.
%
% Uses the same rule as the clean pipeline: fit a linear model from
% observed CSD correlation to the requested information-flow outcome, then
% choose the CSD correlation where the lower 95% confidence bound crosses 0.

if nargin < 2 || isempty(outcomeName)
    outcomeName = "info_observed_corr";
end
if nargin < 3 || isempty(fallbackThreshold)
    fallbackThreshold = 0.30;
end
outcomeName = string(outcomeName);

assert(istable(T), 'Input T must be a validation summary table.');
assert(ismember("csd_observed_abs_corr", string(T.Properties.VariableNames)), ...
    'Validation table must contain csd_observed_abs_corr.');
assert(ismember(outcomeName, string(T.Properties.VariableNames)), ...
    'Threshold outcome %s is not in validation summary table.', char(outcomeName));

x = T.csd_observed_abs_corr;
y = T.(char(outcomeName));
good = isfinite(x) & isfinite(y);
if nnz(good) < 3 || range(x(good)) == 0
    thresholdFit.threshold = fallbackThreshold;
    thresholdFit.status = "fallback_too_few_points";
    thresholdFit.fallback_threshold = fallbackThreshold;
    thresholdFit.outcome = outcomeName;
    thresholdFit.n_points = nnz(good);
    thresholdFit.fit_table = table();
    return;
end

mdl = fitlm(x(good), y(good));
xGrid = linspace(max(0, min(x(good))), min(1, max(x(good))), 400).';
[yHat, yCI] = predict(mdl, xGrid, 'Alpha', 0.05, 'Prediction', 'curve');
lower95 = yCI(:,1);
idx = find(lower95 >= 0, 1, 'first');

if isempty(idx)
    threshold = fallbackThreshold;
    status = "fallback_lower95_never_crosses_zero";
elseif idx == 1
    threshold = xGrid(1);
    status = "derived_lower95_already_nonnegative";
else
    x0 = xGrid(idx-1);
    x1 = xGrid(idx);
    y0 = lower95(idx-1);
    y1 = lower95(idx);
    threshold = x0 + (0 - y0) * (x1 - x0) / max(y1 - y0, eps);
    status = "derived_lower95_crosses_zero";
end
threshold = max(0, min(1, threshold));

thresholdFit.threshold = threshold;
thresholdFit.status = status;
thresholdFit.fallback_threshold = fallbackThreshold;
thresholdFit.outcome = outcomeName;
thresholdFit.n_points = nnz(good);
thresholdFit.model = mdl;
thresholdFit.fit_table = table(xGrid, yHat, lower95, yCI(:,2), ...
    'VariableNames', {'csd_observed_abs_corr', 'predicted_infoflow_corr', ...
    'lower95', 'upper95'});
end
