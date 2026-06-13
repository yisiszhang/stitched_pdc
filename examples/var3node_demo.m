% var3node_demo.m
% Synthetic VAR(1) demo for stitched nonparametric PDC.
%
% Network truth:
%   X3 drives X1 and X2, but there is no direct X1 <-> X2 edge.
%
% The point of the demo is to separate two failure modes:
%   1. Pairwise-only PDC on a hidden-common-driver pair can create a
%      spurious X1-X2 connection. Stitched 3-node PDC suppresses this when
%      all three pairwise CSD blocks are available.
%   2. If only X1-X2 and X2-X3 are observed, NNM completion can still infer
%      a 3-node PDC, while naive zero-fill completion is a useful bad control.

clear; close all; clc;
setup('compile');

rng(7);

% -------------------------------------------------------------------------
% VAR(1): columns are sources at t-1, rows are targets at t.
% -------------------------------------------------------------------------
a11 = 0.30;  % X1 self-coupling
a22 = 0.30;  % X2 self-coupling
a33 = 0.80;  % persistent hidden/common driver
a13 = 0.60;  % X3 -> X1
a23 = 0.60;  % X3 -> X2

%       X1    X2    X3
A = [a11,   0,    a13;   % X1
       0,   a22,  a23;   % X2
       0,   0,    a33];  % X3

C = eye(3) * 0.1;
w = zeros(3, 1);

N = 3;
n = 10000;
ndisc = 1000;
labels = {'X1','X2','X3'};

recFull = {1:N};
recAllPairs = {[1 2], [2 3], [1 3]};
recChain = {[1 2], [2 3]};

% -------------------------------------------------------------------------
% Smoother spectral settings for publication-style curves.
% -------------------------------------------------------------------------
params.fs = 1;
params.win = bartlett(16);
params.nov = 10;
params.nfft = 256;
params.lambda = 0;
params.regularizer = 'eigfloor';

% Full simultaneous reference.
xFull = {arsim(w, A, C, n, ndisc)};

% Independent pair recordings.
xAllPairs = local_simulate_blocks(w, A, C, n, ndisc, recAllPairs);
xChain = xAllPairs(1:2);

% -------------------------------------------------------------------------
% Estimate spectra and PDC.
% -------------------------------------------------------------------------
params.method = 'none';
[~, Sfull, SregFull, f] = reconstruct_inversepsd(xFull, recFull, params);
pdcFull = nonparam_pdc_H(SregFull, f);

params.method = 'nnm';
[~, SallPairs, SregAllPairs, fAllPairs] = reconstruct_inversepsd(xAllPairs, recAllPairs, params);
pdcAllPairs = nonparam_pdc_H(SregAllPairs, fAllPairs);

[~, Schain, SregChain, fChain] = reconstruct_inversepsd(xChain, recChain, params);
pdcChain = nonparam_pdc_H(SregChain, fChain);

params.method = 'naive';
[~, SnaiveChain, SregNaiveChain, fNaiveChain] = reconstruct_inversepsd(xChain, recChain, params);
pdcNaiveChain = nonparam_pdc_H(SregNaiveChain, fNaiveChain);

% Pairwise controls: estimate PDC inside each observed pair only, then embed
% those curves back into a 3-node display. These controls intentionally
% ignore the third node and therefore expose hidden-common-driver artifacts.
pdcPairAll = local_pairwise_pdc_controls(xAllPairs, recAllPairs, params, N);
pdcPairChain = local_pairwise_pdc_controls(xChain, recChain, params, N);

pdcTruth = local_var1_pdc(A, f, params.fs);

% -------------------------------------------------------------------------
% Figure 1: all three pairwise blocks are available.
% -------------------------------------------------------------------------
fig1 = figure('Color', 'w', 'Position', [60 80 980 780]);
curves1 = { ...
    struct('name', 'Truth',              'pdc', pdcTruth,    'style', '-',  'width', 2.0, 'color', [0.05 0.05 0.05]), ...
    struct('name', 'Full simultaneous',  'pdc', pdcFull,     'style', '-',  'width', 1.5, 'color', [0.00 0.42 0.70]), ...
    struct('name', 'Stitched all pairs', 'pdc', pdcAllPairs, 'style', '--', 'width', 1.8, 'color', [0.84 0.37 0.00]), ...
    struct('name', 'Pairwise-only ctrl', 'pdc', pdcPairAll,  'style', ':',  'width', 2.0, 'color', [0.55 0.16 0.55])};
local_plot_pdc_grid(f, curves1, labels, ...
    'All pairwise CSD blocks: stitching mitigates the spurious X1-X2 edge');

% -------------------------------------------------------------------------
% Figure 2: only chain blocks are observed: X1-X2 and X2-X3.
% -------------------------------------------------------------------------
fig2 = figure('Color', 'w', 'Position', [90 110 980 780]);
curves2 = { ...
    struct('name', 'Truth',                'pdc', pdcTruth,       'style', '-',  'width', 2.0, 'color', [0.05 0.05 0.05]), ...
    struct('name', 'Full simultaneous',    'pdc', pdcFull,        'style', '-',  'width', 1.5, 'color', [0.00 0.42 0.70]), ...
    struct('name', 'Stitched chain, NNM',  'pdc', pdcChain,       'style', '--', 'width', 1.8, 'color', [0.84 0.37 0.00]), ...
    struct('name', 'Naive zero-fill ctrl', 'pdc', pdcNaiveChain,  'style', '-.', 'width', 1.4, 'color', [0.62 0.45 0.12]), ...
    struct('name', 'Pairwise-only ctrl',   'pdc', pdcPairChain,   'style', ':',  'width', 2.0, 'color', [0.55 0.16 0.55])};
local_plot_pdc_grid(f, curves2, labels, ...
    'Chain observations only: X1-X3 is imputed before 3-node PDC');

% -------------------------------------------------------------------------
% Numeric summary for the critical non-edge X1 <-> X2.
% -------------------------------------------------------------------------
fprintf('\nMean |PDC| over %.2f-%.2f cycles/sample for the X1-X2 non-edge\n', ...
    min(f), max(f));
local_print_nonedge_summary('Truth', pdcTruth);
local_print_nonedge_summary('Full simultaneous', pdcFull);
local_print_nonedge_summary('Stitched all pairs', pdcAllPairs);
local_print_nonedge_summary('Pairwise all-pairs control', pdcPairAll);
local_print_nonedge_summary('Stitched chain, NNM', pdcChain);
local_print_nonedge_summary('Naive chain control', pdcNaiveChain);
local_print_nonedge_summary('Pairwise chain control', pdcPairChain);

% Keep handles in the workspace for interactive inspection.
demo = struct();
demo.A = A;
demo.f = f;
demo.Sfull = Sfull;
demo.SallPairs = SallPairs;
demo.Schain = Schain;
demo.SnaiveChain = SnaiveChain;
demo.pdcTruth = pdcTruth;
demo.pdcFull = pdcFull;
demo.pdcAllPairs = pdcAllPairs;
demo.pdcChain = pdcChain;
demo.pdcNaiveChain = pdcNaiveChain;
demo.pdcPairAll = pdcPairAll;
demo.pdcPairChain = pdcPairChain;
demo.fig1 = fig1;
demo.fig2 = fig2;

function xBlocks = local_simulate_blocks(w, A, C, n, ndisc, recset)
xBlocks = cell(1, numel(recset));
for u = 1:numel(recset)
    v = arsim(w, A, C, n, ndisc);
    xBlocks{u} = v(:, recset{u});
end
end

function pdcPair = local_pairwise_pdc_controls(xBlocks, recset, params, N)
pdcPair = [];
for u = 1:numel(recset)
    localParams = params;
    localParams.method = 'none';
    [~, ~, SregPair, fPair] = reconstruct_inversepsd({xBlocks{u}}, {1:numel(recset{u})}, localParams);
    pdcLocal = nonparam_pdc_H(SregPair, fPair);

    if isempty(pdcPair)
        pdcPair = nan(N, N, numel(fPair));
    end

    idx = recset{u};
    pdcPair(idx, idx, :) = pdcLocal;
end
end

function pdc = local_var1_pdc(A, f, fs)
N = size(A, 1);
nf = numel(f);
pdc = zeros(N, N, nf);
for k = 1:nf
    Af = eye(N) - A * exp(-1i * 2 * pi * f(k) / fs);
    denom = sqrt(sum(abs(Af).^2, 1));
    pdc(:,:,k) = Af ./ repmat(denom, N, 1);
end
end

function local_plot_pdc_grid(f, curves, labels, titleText)
N = numel(labels);
lineHandles = gobjects(1, numel(curves));
for i = 1:N
    for j = 1:N
        subplot(N, N, (i-1)*N + j);
        hold on;
        for c = 1:numel(curves)
            y = squeeze(abs(curves{c}.pdc(i,j,:)));
            if all(isnan(y))
                continue;
            end
            h = plot(f, y, curves{c}.style, ...
                'LineWidth', curves{c}.width, ...
                'Color', curves{c}.color);
            if i == 1 && j == 1
                lineHandles(c) = h;
            end
        end
        ylim([0 1]);
        xlim([min(f) max(f)]);
        box off;
        set(gca, 'TickDir', 'out', 'FontName', 'Helvetica', 'FontSize', 9);
        if j == 1
            ylabel(sprintf('to %s', labels{i}));
        end
        if i == N
            xlabel(sprintf('from %s', labels{j}));
        end
        if i == 1
            title(sprintf('from %s', labels{j}));
        end
        if j == 1
            text(min(f), 0.92, sprintf('to %s', labels{i}), ...
                'FontWeight', 'bold', 'FontSize', 9);
        end
    end
end
validLegend = isgraphics(lineHandles);
legendAx = axes('Position', [0.10 0.01 0.80 0.04], 'Visible', 'off');
legend(legendAx, lineHandles(validLegend), ...
    cellfun(@(s) s.name, curves(validLegend), 'UniformOutput', false), ...
    'Location', 'north', 'Orientation', 'horizontal', 'Box', 'off');
sgtitle(titleText, 'FontWeight', 'bold');
end

function local_print_nonedge_summary(name, pdc)
% Direction convention is column -> row, so X2->X1 is (1,2) and X1->X2 is (2,1).
x2_to_x1 = mean(abs(squeeze(pdc(1,2,:))), 'omitnan');
x1_to_x2 = mean(abs(squeeze(pdc(2,1,:))), 'omitnan');
fprintf('  %-28s X2->X1 %.3f   X1->X2 %.3f\n', name, x2_to_x1, x1_to_x2);
end
