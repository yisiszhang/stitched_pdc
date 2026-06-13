%% RUN_IBL_METHOD_VERIFICATION
% Three-step verification pipeline for IBL stitched PDC analysis.
%
% Outputs are written to:
%   ibl_output/method_verification
%
% Steps:
%   1. Within-session split/subset stitching validation.
%   2. Cross-session source/target validation and CSD/PDC relationship.
%   3. Reliable full-dataset stitching summary and Allen tracing comparison.

clearvars;
close all;

repoRoot = pwd;
setup;

cfg = ibl_default_config();
cfg.verbose = true;

outDir = fullfile(cfg.output_root, 'method_verification');
figDir = fullfile(outDir, 'figures');
if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end
if exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

qc = local_load_qc(cfg);

% -------------------------------------------------------------------------
% User-facing parameters. Change these in one place.
% -------------------------------------------------------------------------
band = [1 80];
allenCsv = fullfile(cfg.output_root, 'allen_tracing', ...
    'allen_tracing_matrix__normalized_projection_volume__max.csv');

withinFile = fullfile(cfg.output_root, 'validation', ...
    'within_session_subset_stitching_coherence_pdc_20repeats.mat');
withinShuffleFile = fullfile(cfg.output_root, 'validation', ...
    'within_session_subset_stitching_coherence_pdc_20repeats_shuffle_null.mat');
crossFile = fullfile(cfg.output_root, 'validation', ...
    'cross_session_filtered_overlap_pdc.mat');
crossAnalysisFile = fullfile(cfg.output_root, 'validation', ...
    'cross_session_filtered_overlap_analysis.mat');
planFile = fullfile(cfg.output_root, 'reliable_session_component_plan.mat');
resultFile = fullfile(cfg.output_root, 'reliable_component_result', ...
    'stitched_pdc_reliable_component.mat');
allenCmpFile = fullfile(cfg.output_root, 'reliable_component_result', ...
    'allen_comparison_reliable_component.mat');

% Step 3 stitching behavior. When true, the script re-runs stitching using
% exactly the sessions selected by the CSD-reproducible reliable component
% plan. This avoids accidentally comparing an older stitched result.
rerunFinalStitching = true;
finalResultDir = fullfile(cfg.output_root, 'reliable_component_result');
finalResultTag = 'reliable_component_rerun';
finalStitchFile = fullfile(finalResultDir, ['stitched_pdc_' finalResultTag '.mat']);
finalAllenCmpFile = fullfile(finalResultDir, ['allen_comparison_' finalResultTag '.mat']);

% Edit these if you want to change final stitching settings.
finalStitchNormalize = 'coherence';    % 'coherence' or 'none'
finalStitchRegularizer = 'eigfloor';   % 'eigfloor' or 'glasso'
finalStitchLambda = 0;
finalStitchParallel = false;

fprintf('\n[verification] Step 1: within-session validation\n');
within = local_step1_within(withinFile, withinShuffleFile, outDir, figDir);

fprintf('\n[verification] Step 2: cross-session validation\n');
cross = local_step2_cross(crossFile, crossAnalysisFile, outDir, figDir);

fprintf('\n[verification] Step 3: final reliable stitching + Allen comparison\n');
final = local_step3_final(cfg, qc, resultFile, planFile, allenCmpFile, allenCsv, band, ...
    outDir, figDir, rerunFinalStitching, finalStitchFile, finalAllenCmpFile, ...
    finalStitchNormalize, finalStitchRegularizer, finalStitchLambda, finalStitchParallel);

verification = struct();
verification.within = within;
verification.cross = cross;
verification.final = final;
verification.band = band;
verification.created_at = string(datetime('now'));
save(fullfile(outDir, 'ibl_method_verification_summary.mat'), ...
    'verification', 'cfg', '-v7.3');

fprintf('\n[verification] Complete. Outputs saved to:\n  %s\n', outDir);


%% ------------------------------------------------------------------------
% Step 1
% -------------------------------------------------------------------------
function out = local_step1_within(withinFile, shuffleFile, outDir, figDir)
assert(exist(withinFile, 'file') == 2, ...
    'Within-session validation file not found: %s', withinFile);
tmp = load(withinFile, 'val');
val = tmp.val;
T = val.summary_table;

records = val.records;
[T.info_all_mse, T.info_observed_mse, T.info_completed_mse] = ...
    local_infoflow_mse_from_records(records);
[T.pdc_all_mse, T.pdc_observed_mse, T.pdc_completed_mse, pdcNote] = ...
    local_pdc_mse_from_records(records);
[T.pc1_similarity_median, pcNote] = local_pc1_similarity_from_records(records);

summary = local_group_summary(T, "within_session");
writetable(T, fullfile(outDir, 'step1_within_session_repeat_metrics.csv'));
writetable(summary, fullfile(outDir, 'step1_within_session_summary.csv'));

shuffle = [];
if exist(shuffleFile, 'file') == 2
    tmp = load(shuffleFile, 'shuffle');
    shuffle = tmp.shuffle;
    if isfield(shuffle, 'summary_table')
        writetable(shuffle.summary_table, fullfile(outDir, 'step1_within_session_shuffle_summary.csv'));
    end
end

fig = local_new_figure([7.2 3.0]);
tiledlayout(1,3, 'TileSpacing', 'compact', 'Padding', 'compact');
local_box_panel(["All","Observed","Unobserved"], ...
    {T.info_all_corr, T.info_observed_corr, T.info_completed_corr}, ...
    'Information-flow correlation', 'a');
local_box_panel(["All","Observed","Unobserved"], ...
    {T.info_all_mse, T.info_observed_mse, T.info_completed_mse}, ...
    'Information-flow MSE', 'b');
local_box_panel(["All","Observed","Unobserved"], ...
    {T.csd_all_abs_corr, T.csd_observed_abs_corr, T.csd_completed_abs_corr}, ...
    'CSD |corr|', 'c');
exportgraphics(fig, fullfile(figDir, 'fig_step1_within_session_validation.pdf'), 'ContentType', 'vector');
exportgraphics(fig, fullfile(figDir, 'fig_step1_within_session_validation.png'), 'Resolution', 450);

out.file = withinFile;
out.table = T;
out.summary = summary;
out.shuffle = shuffle;
out.pdc_mse_note = pdcNote;
out.pc1_similarity_note = pcNote;
save(fullfile(outDir, 'step1_within_session_validation.mat'), 'out', '-v7.3');
end


%% ------------------------------------------------------------------------
% Step 2
% -------------------------------------------------------------------------
function out = local_step2_cross(crossFile, crossAnalysisFile, outDir, figDir)
assert(exist(crossFile, 'file') == 2, ...
    'Cross-session validation file not found: %s', crossFile);
tmp = load(crossFile, 'val');
val = tmp.val;
T = val.summary_table;

if exist(crossAnalysisFile, 'file') == 2
    tmp = load(crossAnalysisFile);
    if isfield(tmp, 'out')
        ana = tmp.out;
    elseif isfield(tmp, 'ana')
        ana = tmp.ana;
    else
        ana = ibl_analyze_cross_session_validation(val, 'MakeFigure', false);
    end
else
    ana = ibl_analyze_cross_session_validation(val, 'MakeFigure', false);
end

pcTable = local_cross_pc_similarity_table(val);
T = local_join_cross_pc_summary(T, pcTable);

writetable(T, fullfile(outDir, 'step2_cross_session_validation_metrics.csv'));
if isfield(ana, 'relationship_table')
    writetable(ana.relationship_table, fullfile(outDir, 'step2_cross_session_relationship_table.csv'));
end
if isfield(ana, 'selected_table')
    writetable(ana.selected_table, fullfile(outDir, 'step2_cross_session_selected_candidates.csv'));
end
if ~isempty(pcTable)
    writetable(pcTable, fullfile(outDir, 'step2_cross_session_pc_similarity.csv'));
end

fig = local_new_figure([7.2 3.0]);
tiledlayout(1,3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
local_scatter_with_fit(T.csd_observed_abs_corr, T.info_observed_corr, ...
    'Observed CSD |corr|', 'Observed info-flow corr', 'a');

nexttile;
local_scatter_with_fit(T.csd_observed_abs_corr, T.info_all_corr, ...
    'Observed CSD |corr|', 'All-pair info-flow corr', 'b');

nexttile;
if any(isfinite(T.pc1_feature_corr_median))
    local_scatter_with_fit(T.pc1_feature_corr_median, T.csd_observed_abs_corr, ...
        'Median PC1/latent similarity', 'Observed CSD |corr|', 'c');
else
    axis off;
    text(0.05, 0.5, 'PC1 similarity not available in saved validation records.', ...
        'FontName', 'Arial', 'FontSize', 8);
    local_panel_label('c');
end
exportgraphics(fig, fullfile(figDir, 'fig_step2_cross_session_validation.pdf'), 'ContentType', 'vector');
exportgraphics(fig, fullfile(figDir, 'fig_step2_cross_session_validation.png'), 'Resolution', 450);

out.file = crossFile;
out.table = T;
out.analysis = ana;
out.pc_similarity_table = pcTable;
save(fullfile(outDir, 'step2_cross_session_validation.mat'), 'out', '-v7.3');
end


%% ------------------------------------------------------------------------
% Step 3
% -------------------------------------------------------------------------
function out = local_step3_final(cfg, qc, resultFile, planFile, allenCmpFile, allenCsv, band, ...
    outDir, figDir, rerunFinalStitching, finalStitchFile, finalAllenCmpFile, ...
    finalStitchNormalize, finalStitchRegularizer, finalStitchLambda, finalStitchParallel)
assert(exist(planFile, 'file') == 2, ...
    'Reliable session plan file not found: %s', planFile);
tmpPlan = load(planFile, 'plan');
plan = tmpPlan.plan;

if rerunFinalStitching
    if exist(fileparts(finalStitchFile), 'dir') ~= 7
        mkdir(fileparts(finalStitchFile));
    end
    cfgFinal = cfg;
    cfgFinal.session_include_ids = string(plan.selected_session_ids(:));
    cfgFinal.stitch_file = finalStitchFile;
    cfgFinal.stitch.normalize = finalStitchNormalize;
    cfgFinal.stitch.regularizer = finalStitchRegularizer;
    cfgFinal.stitch.lambda = finalStitchLambda;
    cfgFinal.stitch.parallel = finalStitchParallel;

    fprintf('[verification] Re-running final stitching with %d CSD-reliable sessions\n', ...
        numel(cfgFinal.session_include_ids));
    fprintf('[verification] normalize=%s regularizer=%s lambda=%g parallel=%d\n', ...
        cfgFinal.stitch.normalize, cfgFinal.stitch.regularizer, ...
        cfgFinal.stitch.lambda, cfgFinal.stitch.parallel);

    result = ibl_stitch_saved_spectra(cfgFinal, qc);
    save(finalStitchFile, 'result', 'plan', 'cfgFinal', '-v7.3');
    resultFile = finalStitchFile;
    allenCmpFile = finalAllenCmpFile;
else
    assert(exist(resultFile, 'file') == 2, ...
        'Reliable stitched result file not found: %s', resultFile);
    tmp = load(resultFile, 'result');
    result = tmp.result;
end

fMask = result.freqs >= band(1) & result.freqs <= min(band(2), result.freqs(end));
infoflow = ibl_pdc_to_infoflow(result.PDC(:,:,fMask));
observationMask = local_observation_mask(result);

if ~rerunFinalStitching && exist(allenCmpFile, 'file') == 2
    tmpCmp = load(allenCmpFile, 'cmp');
    cmp = tmpCmp.cmp;
else
    assert(exist(allenCsv, 'file') == 2, 'Allen tracing CSV not found: %s', allenCsv);
    cmp = ibl_compare_allen_tracing(result, allenCsv, ...
        'Band', band, ...
        'Order', 'cluster', ...
        'RestrictToAllenSupported', true, ...
        'AllenPositiveThreshold', 0, ...
        'ObservationMask', observationMask, ...
        'MakeMatrixFigure', false, ...
        'MakeRocFigure', false);
    save(allenCmpFile, 'cmp', 'result', 'plan', 'band', '-v7.3');
end

areaSupport = table(string(result.area_names(:)), local_area_support(result), ...
    'VariableNames', {'area', 'n_sessions'});
writetable(areaSupport, fullfile(outDir, 'step3_final_area_support.csv'));

cmpSummary = local_cmp_summary_table(cmp);
writetable(cmpSummary, fullfile(outDir, 'step3_allen_comparison_summary.csv'));

Tinfo = array2table(infoflow, ...
    'VariableNames', matlab.lang.makeValidName(cellstr(result.area_names)));
Tinfo = addvars(Tinfo, string(result.area_names(:)), 'Before', 1, 'NewVariableNames', 'target_area');
writetable(Tinfo, fullfile(outDir, 'step3_infoflow_matrix.csv'));

fig = local_new_figure([7.2 6.0]);
tiledlayout(2,2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
local_matrix_plot(infoflow, string(result.area_names(:)), ...
    sprintf('Information flow %g-%g Hz', band(1), band(2)), 'a');

nexttile;
local_matrix_plot(double(observationMask), string(result.area_names(:)), ...
    'Direct co-observation mask', 'b');

nexttile;
if isfield(cmp, 'infoflow_plot') && isfield(cmp, 'allen_plot')
    local_matrix_plot(cmp.infoflow_plot, string(cmp.shared_areas_plot(:)), ...
        'Allen-supported info flow', 'c');
else
    axis off; local_panel_label('c');
end

nexttile;
if isfield(cmp, 'roc') && isfield(cmp.roc, 'fpr')
    plot(cmp.roc.fpr, cmp.roc.tpr, 'Color', [0.05 0.05 0.05], 'LineWidth', 1.6); hold on;
    plot([0 1], [0 1], '--', 'Color', [0.65 0.65 0.65], 'LineWidth', 0.8);
    axis square; box off; grid on;
    xlabel('False positive rate');
    ylabel('True positive rate');
    title(sprintf('Allen ROC-AUC = %.2f', cmp.roc_auc));
    local_panel_label('d');
else
    axis off; local_panel_label('d');
end
exportgraphics(fig, fullfile(figDir, 'fig_step3_final_stitching_allen.pdf'), 'ContentType', 'vector');
exportgraphics(fig, fullfile(figDir, 'fig_step3_final_stitching_allen.png'), 'Resolution', 450);

out.result_file = resultFile;
out.plan_file = planFile;
out.allen_comparison_file = allenCmpFile;
out.cmp = cmp;
out.cmp_summary = cmpSummary;
out.infoflow = infoflow;
out.observation_mask = observationMask;
out.plan = plan;
save(fullfile(outDir, 'step3_final_stitching_allen.mat'), 'out', '-v7.3');
end


%% ------------------------------------------------------------------------
% Metric helpers
% -------------------------------------------------------------------------
function [mAll, mObs, mComp] = local_infoflow_mse_from_records(records)
n = numel(records);
mAll = nan(n,1); mObs = nan(n,1); mComp = nan(n,1);
for i = 1:n
    if ~isfield(records(i), 'infoflow_stitch') || ~isfield(records(i), 'infoflow_full_control')
        continue;
    end
    A = records(i).infoflow_stitch;
    B = records(i).infoflow_full_control;
    obs = ~records(i).missing_mask;
    off = ~eye(size(obs));
    observed = obs & off;
    completed = ~obs & off;
    mAll(i) = local_mse(A(off), B(off));
    mObs(i) = local_mse(A(observed), B(observed));
    mComp(i) = local_mse(A(completed), B(completed));
end
end


function [mAll, mObs, mComp, note] = local_pdc_mse_from_records(records)
n = numel(records);
mAll = nan(n,1); mObs = nan(n,1); mComp = nan(n,1);
note = "PDC arrays were not stored in the saved validation records; PDC off-diagonal MSE is unavailable. Use info-flow MSE as saved-data proxy or rerun validation with PDC arrays saved.";
for i = 1:n
    if ~isfield(records(i), 'PDC_stitch') || ~isfield(records(i), 'PDC_full_control')
        continue;
    end
    A = records(i).PDC_stitch;
    B = records(i).PDC_full_control;
    obs = ~records(i).missing_mask;
    off = ~eye(size(obs));
    observed = obs & off;
    completed = ~obs & off;
    mAll(i) = local_pdc_mask_mse(A, B, off);
    mObs(i) = local_pdc_mask_mse(A, B, observed);
    mComp(i) = local_pdc_mask_mse(A, B, completed);
    note = "PDC off-diagonal MSE computed from saved PDC arrays.";
end
end


function [pcSim, note] = local_pc1_similarity_from_records(records)
pcSim = nan(numel(records), 1);
note = "PC1 loading vectors were not stored in the saved within-session records; PC1 similarity is unavailable unless validation is rerun with loadings saved.";
for i = 1:numel(records)
    if isfield(records(i), 'pc1_similarity_median')
        pcSim(i) = records(i).pc1_similarity_median;
        note = "PC1 similarity loaded from saved records.";
    end
end
end


function mse = local_pdc_mask_mse(A, B, mask)
idx = find(mask);
if isempty(idx)
    mse = NaN;
    return;
end
valsA = [];
valsB = [];
for k = 1:numel(idx)
    [i,j] = ind2sub(size(mask), idx(k));
    valsA = [valsA; squeeze(A(i,j,:))]; %#ok<AGROW>
    valsB = [valsB; squeeze(B(i,j,:))]; %#ok<AGROW>
end
mse = local_mse(valsA, valsB);
end


function y = local_mse(a, b)
a = real(a(:)); b = real(b(:));
good = isfinite(a) & isfinite(b);
if ~any(good)
    y = NaN;
else
    d = a(good) - b(good);
    y = mean(d.^2);
end
end


function T = local_group_summary(Tin, label)
vars = ["info_all_corr","info_observed_corr","info_completed_corr", ...
    "info_all_mse","info_observed_mse","info_completed_mse", ...
    "csd_all_abs_corr","csd_observed_abs_corr","csd_completed_abs_corr"];
metric = strings(numel(vars),1);
median_value = nan(numel(vars),1);
iqr_low = nan(numel(vars),1);
iqr_high = nan(numel(vars),1);
n = zeros(numel(vars),1);
for i = 1:numel(vars)
    metric(i) = vars(i);
    if ismember(vars(i), string(Tin.Properties.VariableNames))
        x = Tin.(vars(i));
        x = x(isfinite(x));
        n(i) = numel(x);
        if ~isempty(x)
            median_value(i) = median(x);
            iqr_low(i) = prctile(x, 25);
            iqr_high(i) = prctile(x, 75);
        end
    end
end
stage = repmat(string(label), numel(vars), 1);
T = table(stage, metric, n, median_value, iqr_low, iqr_high);
end


function pcTable = local_cross_pc_similarity_table(val)
records = val.records;
rows = struct('target_session_id', {}, 'source_session_id', {}, ...
    'area', {}, 'pc1_feature_corr', {});
% Saved cross-session records generally do not contain PCA feature vectors.
% If future records include a fingerprint table, this hook will export it.
for i = 1:numel(records)
    if isfield(records(i), 'pc_similarity_table')
        S = records(i).pc_similarity_table;
        for r = 1:height(S)
            row.target_session_id = string(records(i).target_session_id);
            row.source_session_id = string(S.source_session_id(r));
            row.area = string(S.area(r));
            if ismember('feature_corr', string(S.Properties.VariableNames))
                row.pc1_feature_corr = S.feature_corr(r);
            else
                row.pc1_feature_corr = NaN;
            end
            rows(end+1) = row; %#ok<AGROW>
        end
    end
end
if isempty(rows)
    pcTable = table();
else
    pcTable = struct2table(rows);
end
end


function T = local_join_cross_pc_summary(T, pcTable)
T.pc1_feature_corr_median = nan(height(T),1);
if isempty(pcTable)
    return;
end
for i = 1:height(T)
    idx = pcTable.target_session_id == T.target_session_id(i);
    T.pc1_feature_corr_median(i) = median(pcTable.pc1_feature_corr(idx), 'omitnan');
end
end


function obs = local_observation_mask(result)
n = numel(result.area_names);
if isfield(result, 'meta') && isfield(result.meta, 'missing_mask')
    obs = ~logical(result.meta.missing_mask);
elseif isfield(result, 'recset')
    obs = meacount_mat(result.recset, n) > 0;
else
    obs = true(n);
end
obs(1:n+1:end) = false;
end


function support = local_area_support(result)
n = numel(result.area_names);
support = zeros(n,1);
if isfield(result, 'recset')
    for k = 1:numel(result.recset)
        idx = result.recset{k};
        support(idx) = support(idx) + 1;
    end
end
end


function T = local_cmp_summary_table(cmp)
fields = {'pearson_r','spearman_r','valid_pairs','n_positive','n_negative', ...
    'positive_rate','roc_auc','average_precision','precision_at_npos', ...
    'recall_at_npos','observed_valid_pairs','never_observed_valid_pairs', ...
    'observed_roc_auc','never_observed_roc_auc', ...
    'observed_average_precision','never_observed_average_precision'};
metric = strings(numel(fields),1);
value = nan(numel(fields),1);
for i = 1:numel(fields)
    metric(i) = string(fields{i});
    if isfield(cmp, fields{i})
        v = cmp.(fields{i});
        if isnumeric(v) && isscalar(v)
            value(i) = v;
        end
    end
end
T = table(metric, value);
end


function qc = local_load_qc(cfg)
if exist(cfg.pca_qc_file, 'file') == 2
    tmp = load(cfg.pca_qc_file, 'qc');
    qc = tmp.qc;
else
    qc = [];
end
end


%% ------------------------------------------------------------------------
% Plot helpers
% -------------------------------------------------------------------------
function fig = local_new_figure(sizeInches)
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 sizeInches]);
set(fig, 'Renderer', 'painters');
end


function local_box_panel(labels, values, ylab, panel)
nexttile;
hold on;
colors = [0.10 0.32 0.50; 0.20 0.55 0.45; 0.74 0.36 0.18];
for i = 1:numel(values)
    x = values{i};
    x = x(isfinite(x));
    if isempty(x)
        continue;
    end
    swarmX = i + 0.13 * (rand(size(x)) - 0.5);
    scatter(swarmX, x, 12, 'MarkerFaceColor', colors(i,:), ...
        'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.55);
    q = prctile(x, [25 50 75]);
    plot([i-0.22 i+0.22], [q(2) q(2)], '-', 'Color', [0 0 0], 'LineWidth', 1.2);
    rectangle('Position', [i-0.16 q(1) 0.32 q(3)-q(1)], ...
        'EdgeColor', [0 0 0], 'LineWidth', 0.8);
end
set(gca, 'XTick', 1:numel(labels), 'XTickLabel', labels, 'FontName', 'Arial', 'FontSize', 7);
ylabel(ylab, 'FontName', 'Arial', 'FontSize', 8);
box off; grid on;
local_panel_label(panel);
end


function local_scatter_with_fit(x, y, xlab, ylab, panel)
good = isfinite(x) & isfinite(y);
scatter(x(good), y(good), 30, 'MarkerFaceColor', [0.12 0.35 0.55], ...
    'MarkerEdgeColor', 'w', 'LineWidth', 0.4); hold on;
if nnz(good) >= 3
    mdl = fitlm(x(good), y(good));
    xx = linspace(min(x(good)), max(x(good)), 100).';
    [yy, ci] = predict(mdl, xx, 'Prediction', 'curve');
    plot(xx, yy, 'k-', 'LineWidth', 1.2);
    plot(xx, ci(:,1), 'k--', 'LineWidth', 0.8);
    plot(xx, ci(:,2), 'k--', 'LineWidth', 0.8);
    [rho, p] = corr(x(good), y(good), 'type', 'Spearman');
    title(sprintf('\\rho=%.2f, p=%.3g', rho, p), 'FontWeight', 'normal');
end
yline(0, '-', 'Color', [0.75 0.2 0.2], 'LineWidth', 0.8);
xlabel(xlab, 'FontName', 'Arial', 'FontSize', 8);
ylabel(ylab, 'FontName', 'Arial', 'FontSize', 8);
set(gca, 'FontName', 'Arial', 'FontSize', 7);
box off; grid on;
local_panel_label(panel);
end


function local_matrix_plot(M, labels, titleText, panel)
imagesc(M);
axis image;
colormap(gca, local_viridis_colormap(256));
colorbar;
set(gca, 'XTick', 1:numel(labels), 'XTickLabel', labels, ...
    'YTick', 1:numel(labels), 'YTickLabel', labels, ...
    'XTickLabelRotation', 90, 'TickLabelInterpreter', 'none', ...
    'FontName', 'Arial', 'FontSize', 5);
title(titleText, 'FontName', 'Arial', 'FontSize', 8, 'FontWeight', 'normal');
xlabel('Source area', 'FontName', 'Arial', 'FontSize', 8);
ylabel('Target area', 'FontName', 'Arial', 'FontSize', 8);
local_panel_label(panel);
end


function local_panel_label(panel)
text(-0.12, 1.08, panel, 'Units', 'normalized', ...
    'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 10);
end


function cmap = local_nature_colormap()
x = linspace(0, 1, 256).';
c1 = [0.05 0.10 0.25];
c2 = [0.93 0.94 0.88];
c3 = [0.65 0.12 0.08];
cmap = zeros(256,3);
for i = 1:256
    if x(i) < 0.5
        t = x(i) / 0.5;
        cmap(i,:) = (1-t)*c1 + t*c2;
    else
        t = (x(i)-0.5) / 0.5;
        cmap(i,:) = (1-t)*c2 + t*c3;
    end
end
end


function cmap = local_viridis_colormap(n)
if nargin < 1 || isempty(n)
    n = 256;
end
anchors = [ ...
    0.267004 0.004874 0.329415
    0.282623 0.140926 0.457517
    0.253935 0.265254 0.529983
    0.206756 0.371758 0.553117
    0.163625 0.471133 0.558148
    0.127568 0.566949 0.550556
    0.134692 0.658636 0.517649
    0.266941 0.748751 0.440573
    0.477504 0.821444 0.318195
    0.741388 0.873449 0.149561
    0.993248 0.906157 0.143936];
x = linspace(0, 1, size(anchors, 1));
xi = linspace(0, 1, n);
cmap = interp1(x, anchors, xi, 'pchip');
cmap = max(0, min(1, cmap));
end
