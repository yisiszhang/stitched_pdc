%% PLOT_IBL_METHOD_FIGURES_1_2_FROM_CLEAN_RUN
% Recreate method Figures 1 and 2 from the clean-run outputs, matching the
% style of plot_ibl_method_figures_from_saved.m but without using old
% ibl_output/validation paths.

clearvars;
close all;

repoRoot = pwd;
setup;

runRoot = fullfile(repoRoot, 'ibl_output', 'clean_runs', ...
    'clean_maxneurons100_20260515_025635');

outDir = fullfile(runRoot, 'method_figures_from_clean_run');
if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end

withinFile = fullfile(runRoot, '02_within_session_validation', ...
    'within_session_subset_stitching_20repeats.mat');
crossFile = fullfile(runRoot, '05_reliable_plan', ...
    'threshold_cross_session_validation.mat');

assert(exist(withinFile, 'file') == 2, ...
    'Within-session validation file not found: %s', withinFile);
assert(exist(crossFile, 'file') == 2, ...
    'Cross-session validation file not found: %s', crossFile);

fprintf('[clean plot] loading\n  %s\n  %s\n', withinFile, crossFile);
valWithin = local_load_val(withinFile);
valCross = local_load_val(crossFile);

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

local_export(fig1, outDir, 'fig1_within_session_validation_clean_run');

% -------------------------------------------------------------------------
% Figure 2: cross-session CSD reproducibility vs PDC information flow.
% -------------------------------------------------------------------------
Tcross = valCross.summary_table;
fig2 = local_new_figure([3.45 3.1]);
local_scatter_with_fit(Tcross.csd_observed_abs_corr, Tcross.info_observed_corr, ...
    'Observed CSD |corr|', 'Observed info-flow corr', 'a');
local_export(fig2, outDir, 'fig2_cross_session_csd_vs_infoflow_clean_run');

% Source tables.
writetable(Twithin, fullfile(outDir, 'fig1_within_session_source_data_clean_run.csv'));
writetable(Tcross, fullfile(outDir, 'fig2_cross_session_source_data_clean_run.csv'));

fprintf('[clean plot] figures saved to:\n  %s\n', outDir);

function val = local_load_val(filename)
S = load(filename);
if isfield(S, 'val')
    val = S.val;
elseif isfield(S, 'valWithin')
    val = S.valWithin;
elseif isfield(S, 'valThreshold')
    val = S.valThreshold;
else
    names = fieldnames(S);
    val = [];
    for i = 1:numel(names)
        x = S.(names{i});
        if isstruct(x) && isfield(x, 'summary_table')
            val = x;
            break;
        end
    end
    assert(~isempty(val), 'No validation struct with summary_table found in %s.', filename);
end
assert(isfield(val, 'summary_table'), ...
    'Validation file %s does not contain summary_table.', filename);
end

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

function local_panel_label(panel)
text(-0.14, 1.08, panel, 'Units', 'normalized', ...
    'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 10);
end
