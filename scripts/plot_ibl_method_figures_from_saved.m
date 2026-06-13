%% PLOT_IBL_METHOD_FIGURES_FROM_SAVED
% Plot publication-style method verification figures from saved outputs.
%
% Run from repo root:
%   cd('/Users/metis/Documents/MATLAB/stitch_causality/icml/snpdc_package/snpdc_icml')
%   run('scripts/plot_ibl_method_figures_from_saved.m')

clearvars;close all;

repoRoot = pwd;
setup;

cfg = ibl_default_config();

outDir = fullfile(cfg.output_root, 'method_verification', 'figures_from_saved');
if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end

% -------------------------------------------------------------------------
% Saved result files. Edit these if you want to plot another rerun.
% -------------------------------------------------------------------------
withinFile = fullfile(cfg.output_root, 'validation', ...
    'within_session_subset_stitching_coherence_pdc_20repeats.mat');
crossFile = fullfile(cfg.output_root, 'validation', ...
    'cross_session_filtered_overlap_pdc.mat');
infoflowFile = fullfile(cfg.output_root, 'reliable_component_result', ...
    'infoflow_reliable_component.mat');
allenCmpFile = fullfile(cfg.output_root, 'reliable_component_result', ...
    'allen_comparison_reliable_component.mat');
stitchedFile = fullfile(cfg.output_root, 'reliable_component_result', ...
    'stitched_pdc_reliable_component.mat');

fprintf('[plot] loading saved results\n');
within = load(withinFile, 'val');
cross = load(crossFile, 'val');
info = load(infoflowFile);
allen = load(allenCmpFile, 'cmp');
stitched = load(stitchedFile, 'result');

valWithin = within.val;
valCross = cross.val;
cmp = allen.cmp;
result = stitched.result;
if isfield(info, 'infoflow')
    infoflow = info.infoflow;
else
    band = [1 80];
    fMask = result.freqs >= band(1) & result.freqs <= min(band(2), result.freqs(end));
    infoflow = ibl_pdc_to_infoflow(result.PDC(:,:,fMask));
end

% -------------------------------------------------------------------------
% Figure 1: within-session partial-observation validation.
% -------------------------------------------------------------------------
Twithin = local_within_table_with_mse(valWithin);
fig1 = local_new_figure([7.1 2.65]);
tiledlayout(1,3, 'TileSpacing', 'compact', 'Padding', 'compact');

local_box_panel(["All","Observed","Unobserved"], ...
    {Twithin.info_all_corr, Twithin.info_observed_corr, Twithin.info_completed_corr}, ...
    'Info-flow correlation', 'a');
yline(0, '-', 'Color', [0.65 0.65 0.65], 'LineWidth', 0.8);

local_box_panel(["All","Observed","Unobserved"], ...
    {Twithin.info_all_mse, Twithin.info_observed_mse, Twithin.info_completed_mse}, ...
    'Info-flow MSE', 'b');

local_box_panel(["All","Observed","Unobserved"], ...
    {Twithin.csd_all_abs_corr, Twithin.csd_observed_abs_corr, Twithin.csd_completed_abs_corr}, ...
    'CSD |corr|', 'c');
yline(0, '-', 'Color', [0.65 0.65 0.65], 'LineWidth', 0.8);

local_export(fig1, outDir, 'fig1_within_session_validation');

% -------------------------------------------------------------------------
% Figure 2: cross-session CSD reproducibility vs PDC information flow.
% -------------------------------------------------------------------------
Tcross = valCross.summary_table;
fig2 = local_new_figure([3.45 3.1]);
local_scatter_with_fit(Tcross.csd_observed_abs_corr, Tcross.info_observed_corr, ...
    'Observed CSD |corr|', 'Observed info-flow corr', 'a');
local_export(fig2, outDir, 'fig2_cross_session_csd_vs_infoflow');

% -------------------------------------------------------------------------
% Figure 3: final reliable-component stitching and Allen comparison.
% -------------------------------------------------------------------------
obsMask = local_observation_mask(result);
fig3 = local_new_figure([7.1 6.2]);
tiledlayout(2,2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
local_matrix_plot(infoflow, string(result.area_names(:)), ...
    'Stitched information flow', 'a', 'flow', 'scale');

nexttile;
if isfield(cmp, 'allen_plot')
    local_matrix_plot(cmp.allen_plot, string(cmp.shared_areas_plot(:)), ...
        'Allen tracing', 'b', 'allen', 'cutoff');
else
    axis off;
    local_panel_label('b');
    text(0.05, 0.5, 'Allen matrix unavailable', 'FontName', 'Arial', 'FontSize', 8);
end

nexttile;
local_matrix_plot(double(obsMask), string(result.area_names(:)), ...
    'Direct co-observation', 'c', 'mask','orig');

nexttile;
local_roc_panel(cmp, 'd');

local_export(fig3, outDir, 'fig3_final_stitching_allen');

% -------------------------------------------------------------------------
% Export concise figure source tables.
% -------------------------------------------------------------------------
writetable(Twithin, fullfile(outDir, 'fig1_within_session_source_data.csv'));
writetable(Tcross, fullfile(outDir, 'fig2_cross_session_source_data.csv'));
local_write_matrix_csv(infoflow, string(result.area_names(:)), ...
    fullfile(outDir, 'fig3_infoflow_matrix_source_data.csv'));
if isfield(cmp, 'allen_plot')
    local_write_matrix_csv(cmp.allen_plot, string(cmp.shared_areas_plot(:)), ...
        fullfile(outDir, 'fig3_allen_matrix_source_data.csv'));
end
local_write_cmp_summary(cmp, fullfile(outDir, 'fig3_allen_metrics_source_data.csv'));

fprintf('[plot] figures saved to:\n  %s\n', outDir);


%% ------------------------------------------------------------------------
% Data helpers
% -------------------------------------------------------------------------
function T = local_within_table_with_mse(val)
T = val.summary_table;
records = val.records;
n = numel(records);
T.info_all_mse = nan(n,1);
T.info_observed_mse = nan(n,1);
T.info_completed_mse = nan(n,1);
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
    T.info_all_mse(i) = local_mse(A(off), B(off));
    T.info_observed_mse(i) = local_mse(A(observed), B(observed));
    T.info_completed_mse(i) = local_mse(A(completed), B(completed));
end
end


function y = local_mse(a, b)
a = real(a(:));
b = real(b(:));
good = isfinite(a) & isfinite(b);
if ~any(good)
    y = NaN;
else
    d = a(good) - b(good);
    y = mean(d.^2);
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


function local_write_matrix_csv(M, labels, filename)
T = array2table(M, 'VariableNames', matlab.lang.makeValidName(cellstr(labels)));
T = addvars(T, labels(:), 'Before', 1, 'NewVariableNames', 'target_area');
writetable(T, filename);
end


function local_write_cmp_summary(cmp, filename)
fields = {'pearson_r','spearman_r','valid_pairs','n_positive','n_negative', ...
    'positive_rate','roc_auc','average_precision','precision_at_npos', ...
    'recall_at_npos','observed_roc_auc','never_observed_roc_auc', ...
    'observed_average_precision','never_observed_average_precision'};
metric = strings(numel(fields),1);
value = nan(numel(fields),1);
for i = 1:numel(fields)
    metric(i) = string(fields{i});
    if isfield(cmp, fields{i}) && isnumeric(cmp.(fields{i})) && isscalar(cmp.(fields{i}))
        value(i) = cmp.(fields{i});
    end
end
writetable(table(metric, value), filename);
end


%% ------------------------------------------------------------------------
% Plot helpers
% -------------------------------------------------------------------------
function fig = local_new_figure(sizeInches)
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 sizeInches]);
set(fig, 'Renderer', 'painters');
end


function local_export(fig, outDir, stem)
set(fig, 'PaperPositionMode', 'auto');
exportgraphics(fig, fullfile(outDir, [stem '.pdf']), 'ContentType', 'vector');
exportgraphics(fig, fullfile(outDir, [stem '.png']), 'Resolution', 600);
savefig(fig, fullfile(outDir, [stem '.fig']));
end


function local_box_panel(labels, values, ylab, panel)
nexttile;
hold on;
colors = [0.08 0.29 0.46; 0.10 0.48 0.40; 0.72 0.31 0.16];
for i = 1:numel(values)
    x = values{i};
    x = x(isfinite(x));
    if isempty(x)
        continue;
    end
    rng(i);
    swarmX = i + 0.16 * (rand(size(x)) - 0.5);
    scatter(swarmX, x, 14, ...
        'MarkerFaceColor', colors(i,:), ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.50);
    q = prctile(x, [25 50 75]);
    rectangle('Position', [i-0.17 q(1) 0.34 max(q(3)-q(1), eps)], ...
        'EdgeColor', [0.05 0.05 0.05], ...
        'LineWidth', 0.8);
    plot([i-0.22 i+0.22], [q(2) q(2)], '-', ...
        'Color', [0.05 0.05 0.05], 'LineWidth', 1.2);
end
set(gca, 'XTick', 1:numel(labels), ...
    'XTickLabel', labels, ...
    'TickDir', 'out', ...
    'FontName', 'Arial', ...
    'FontSize', 7, ...
    'LineWidth', 0.8);
ylabel(ylab, 'FontName', 'Arial', 'FontSize', 8);
box off;
grid on;
local_panel_label(panel);
end


function local_scatter_with_fit(x, y, xlab, ylab, panel)
good = isfinite(x) & isfinite(y);
x = x(good);
y = y(good);
scatter(x, y, 42, ...
    'MarkerFaceColor', [0.08 0.29 0.46], ...
    'MarkerEdgeColor', 'w', ...
    'LineWidth', 0.45, ...
    'MarkerFaceAlpha', 0.82);
hold on;
if numel(x) >= 3
    mdl = fitlm(x, y);
    xx = linspace(min(x), max(x), 200).';
    [yy, ci] = predict(mdl, xx, 'Prediction', 'curve');
    plot(xx, yy, '-', 'Color', [0.05 0.05 0.05], 'LineWidth', 1.5);
    plot(xx, ci(:,1), '--', 'Color', [0.25 0.25 0.25], 'LineWidth', 0.9);
    plot(xx, ci(:,2), '--', 'Color', [0.25 0.25 0.25], 'LineWidth', 0.9);
    [rho, p] = corr(x, y, 'type', 'Spearman');
    text(0.05, 0.92, sprintf('\\rho = %.2f, p = %.2g', rho, p), ...
        'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 8);
end
yline(0, '-', 'Color', [0.70 0.18 0.14], 'LineWidth', 0.9);
xlabel(xlab, 'FontName', 'Arial', 'FontSize', 8);
ylabel(ylab, 'FontName', 'Arial', 'FontSize', 8);
set(gca, 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 7, 'LineWidth', 0.8);
box off;
grid on;
axis square;
local_panel_label(panel);
end


function local_matrix_plot(M, labels, titleText, panel, mode, scale)
switch lower(string(mode))
    case "allen"
        M = M - diag(diag(M));
end
switch lower(string(scale))
    case "scale"
        imagesc(M,[min(M(:)), 1]); %max(M(:))*0.5]);
    case "cutoff"
        imagesc(M,[200, max(M(:))]);
    otherwise
        imagesc(M);
end
axis image;
switch lower(string(mode))
    case "mask"
%         colormap(gca, local_viridis_colormap(2));
        colormap(gca, [0.94 0.94 0.94; 0.08 0.29 0.46]);
        clim([0 1]);
    otherwise
        colormap(gca, local_viridis_colormap(256));
end
colorbar;
set(gca, 'XTick', 1:numel(labels), ...
    'XTickLabel', labels, ...
    'YTick', 1:numel(labels), ...
    'YTickLabel', labels, ...
    'XTickLabelRotation', 90, ...
    'TickLabelInterpreter', 'none', ...
    'TickDir', 'out', ...
    'FontName', 'Arial', ...
    'FontSize', 4.8, ...
    'LineWidth', 0.6);
title(titleText, 'FontName', 'Arial', 'FontSize', 8, 'FontWeight', 'normal');
xlabel('Source area', 'FontName', 'Arial', 'FontSize', 8);
ylabel('Target area', 'FontName', 'Arial', 'FontSize', 8);
local_panel_label(panel);
end


function local_roc_panel(cmp, panel)
if ~isfield(cmp, 'roc') || ~isfield(cmp.roc, 'fpr')
    axis off;
    local_panel_label(panel);
    return;
end
plot(cmp.roc.fpr, cmp.roc.tpr, ...
    'Color', [0.08 0.29 0.46], 'LineWidth', 1.8);
hold on;
legendLabels = strings(0,1);
legendLabels(end+1) = sprintf('All, AUC %.2f', cmp.roc_auc);

if isfield(cmp, 'observed') && isfield(cmp.observed, 'roc') && ...
        isfield(cmp.observed.roc, 'fpr') && ~isempty(cmp.observed.roc.fpr)
    plot(cmp.observed.roc.fpr, cmp.observed.roc.tpr, ...
        'Color', [0.10 0.48 0.40], 'LineWidth', 1.5);
    legendLabels(end+1) = sprintf('Observed, AUC %.2f', cmp.observed_roc_auc);
end

if isfield(cmp, 'never_observed') && isfield(cmp.never_observed, 'roc') && ...
        isfield(cmp.never_observed.roc, 'fpr') && ~isempty(cmp.never_observed.roc.fpr)
    plot(cmp.never_observed.roc.fpr, cmp.never_observed.roc.tpr, ...
        'Color', [0.72 0.31 0.16], 'LineWidth', 1.5);
    legendLabels(end+1) = sprintf('Unobserved, AUC %.2f', cmp.never_observed_roc_auc);
end

plot([0 1], [0 1], '--', 'Color', [0.65 0.65 0.65], 'LineWidth', 0.9);
axis square;
box off;
grid on;
xlabel('False positive rate', 'FontName', 'Arial', 'FontSize', 8);
ylabel('True positive rate', 'FontName', 'Arial', 'FontSize', 8);
title('Allen tracing ROC', ...
    'FontName', 'Arial', 'FontSize', 8, 'FontWeight', 'normal');
legend(legendLabels, 'Location', 'southeast', 'Box', 'off', ...
    'FontName', 'Arial', 'FontSize', 6);
set(gca, 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 7, 'LineWidth', 0.8);
local_panel_label(panel);
end


function local_panel_label(panel)
text(-0.14, 1.08, panel, 'Units', 'normalized', ...
    'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 10);
end


function cmap = local_scientific_colormap()
x = linspace(0, 1, 256).';
c1 = [0.04 0.10 0.23];
c2 = [0.93 0.94 0.88];
c3 = [0.66 0.12 0.08];

cmap = zeros(256,3);
for i = 1:256
    if x(i) < 0.5
        t = x(i) / 0.5;
        cmap(i,:) = (1 - t) * c1 + t * c2;
    else
        t = (x(i) - 0.5) / 0.5;
        cmap(i,:) = (1 - t) * c2 + t * c3;
    end
end
end


function cmap = local_viridis_colormap(n)
% if nargin < 1 || isempty(n)
%     n = 256;
% end
% anchors = [ ...
%     0.267004 0.004874 0.329415
%     0.282623 0.140926 0.457517
%     0.253935 0.265254 0.529983
%     0.206756 0.371758 0.553117
%     0.163625 0.471133 0.558148
%     0.127568 0.566949 0.550556
%     0.134692 0.658636 0.517649
%     0.266941 0.748751 0.440573
%     0.477504 0.821444 0.318195
%     0.741388 0.873449 0.149561
%     0.993248 0.906157 0.143936];
% x = linspace(0, 1, size(anchors, 1));
% xi = linspace(0, 1, n);
% cmap = interp1(x, anchors, xi, 'pchip');
% cmap = max(0, min(1, cmap));
cmap = othercolor('GnBu7');
end
