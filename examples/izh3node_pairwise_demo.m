% izh3node_pairwise_demo.m
% 3-neuron Izhikevich E/I demo with stitched pair observations.
%
% Observation blocks:
%   [1,2], [2,3], [3,1]
%
% Upper bound:
%   full simultaneous observation of all three neurons
%
% Requirements:
%   - Signal Processing Toolbox (for cpsd)
%   - FieldTrip Wilson factorization is bundled in external/

clear; close all; clc;
rng(7);

% setup('compile');

labels = {'N1 (E)', 'N2 (E)', 'N3 (I)'};
recset = {[1 2], [2 3], [3 1]};
N = 3;

% Fixed 2E/1I network. Columns are pre-synaptic neurons; rows are post-synaptic.
W = [0   10  -18; ...
     14   0  -16; ...
      9  12    0];

% Izhikevich parameters: first two excitatory, third inhibitory.
net.a = [0.02; 0.02; 0.06];
net.b = [0.20; 0.20; 0.225];
net.c = [-64.9; -63.7; -65.0];
net.d = [1.0; 1.0; 1.0];
net.W = W;

sim_ms = 200000;
input_mean = 4.0;
noise_amp = 18.0;

% Current Chronux multitaper settings kept for reference:
% chronux.Fs = 1000;
% chronux.tapers = [5 9];
% chronux.pad = 0;
% chronux.fpass = [0 chronux.Fs/2];
% chronux.err = [0 0];
% chronux.trialave = 0;
% chronux.win = 4.096;
% chronux.maxiter = 100;
% chronux.tol = 1e-3;

spec.fs = 1000 / 0.1;
spec.nfft = 2^9;
spec.win = bartlett(2^9);
spec.nov = round(0.875 * numel(spec.win));
spec.fmax = 1000;
spec.maxiter = 100;
spec.tol = 1e-3;

% Full simultaneous observation.
[v_full, dt] = simulate_izhikevich(net, sim_ms, input_mean, noise_amp);
data_full = voltage_to_point_process(v_full, dt);
assert_nonempty_spikes(data_full, 'full observation');
[S_full, f] = cpsd_voltage_spectrum(v_full, spec);
keepF = f <= spec.fmax;
f = f(keepF);
S_full = S_full(:,:,keepF);
pdc_full = nonparam_pdc_H(S_full, f, 'maxiter', spec.maxiter, 'tol', spec.tol);
flow_full = pdc_to_infoflow(pdc_full);

% Stitched pair observations from separate trials.
blockSpectra = cell(1, numel(recset));
for u = 1:numel(recset)
    [v_block, dt_block] = simulate_izhikevich(net, sim_ms, input_mean, noise_amp);
    data_block = voltage_to_point_process(v_block, dt_block);
    assert_nonempty_spikes(data_block, sprintf('stitched block [%d %d]', recset{u}(1), recset{u}(2)));
    [blockSpectra{u}, f_block] = cpsd_voltage_spectrum(v_block(recset{u}, :), spec);
    f_block = f_block(keepF);
    blockSpectra{u} = blockSpectra{u}(:,:,keepF);

    if ~isequal(size(f_block), size(f)) || max(abs(f_block(:) - f(:))) > 1e-12
        error('Frequency grids from cpsd do not match across observation blocks.');
    end
end

S_stitched = stitch_block_spectra(blockSpectra, recset, N);
pdc_stitched = nonparam_pdc_H(S_stitched, f, 'maxiter', spec.maxiter, 'tol', spec.tol);
flow_stitched = pdc_to_infoflow(pdc_stitched);

% Ground-truth flow is directed from column (source) to row (target).
firing_rates = cellfun(@numel, {data_full.times})' / (sim_ms / 1000);
truth_flow = abs(W) .* repmat(firing_rates', N, 1);
truth_flow(1:N+1:end) = 0;

% --------------------------
% Figure 1: PDC comparison
% --------------------------
figure('Position', [60 60 850 760], 'Color', 'w');
for i = 1:N
    for j = 1:N
        subplot(N, N, (i-1)*N + j);
        plot(f, squeeze(abs(pdc_full(i,j,:))), 'LineWidth', 1.5); hold on;
        plot(f, squeeze(abs(pdc_stitched(i,j,:))), '--', 'LineWidth', 1.5);
        xlim([f(1) f(end)]);
        if i == j
            ylim([0 1]);
        else
            ylim([0 0.1]);
        end
        box off;
        if j == 1
            ylabel(labels{i}, 'Interpreter', 'none');
        end
        if i == N
            xlabel(labels{j}, 'Interpreter', 'none');
        end
        title(sprintf('%s -> %s', labels{j}, labels{i}), 'Interpreter', 'none');
    end
end
legend({'Full observation', 'Stitched'}, 'Position', [0.72 0.01 0.22 0.05]);
sgtitle('3-node Izhikevich PDC: full observation vs. stitched pair observations');

% --------------------------------------------
% Figure 2: Information-flow matrix comparison
% --------------------------------------------
flow_full(1:N+1:end) = 0;
flow_stitched(1:N+1:end) = 0;

flow_corr_full = offdiag_corr(truth_flow, flow_full);
flow_corr_stitched = offdiag_corr(truth_flow, flow_stitched);
fprintf('\nCorrelation between synaptic-input truth and information flow (off diagonal):\n');
fprintf('  Full observation: r = %.3f\n', flow_corr_full);
fprintf('  Stitched:         r = %.3f\n', flow_corr_stitched);

figure('Position', [120 120 980 300], 'Color', 'w');

subplot(1,3,1);
plot_flow_matrix(truth_flow, labels);
title('Truth (source = column, target = row)');

subplot(1,3,2);
plot_flow_matrix(flow_full, labels);
title(sprintf('Full observation (r=%.2f)', flow_corr_full));

subplot(1,3,3);
plot_flow_matrix(flow_stitched, labels);
title(sprintf('Stitched (r=%.2f)', flow_corr_stitched));

sgtitle('Information flow (diagonal removed)');


function [v_trace, dt] = simulate_izhikevich(net, T_ms, input_mean, noise_amp)
dt = 0.1;
steps = ceil(T_ms / dt);
N = numel(net.a);

v = -65 * ones(N,1);
u = net.b .* v;
v_trace = zeros(N, steps);

delay_ms = 2;
delay_steps = max(1, round(delay_ms / dt));
spike_buffer = zeros(N, delay_steps);

for t = 1:steps
    I_bg = input_mean + noise_amp * randn(N,1);

    fired = find(v >= 30);
    if ~isempty(fired)
        v(fired) = net.c(fired);
        u(fired) = u(fired) + net.d(fired);
    end

    spike_buffer(:, 1:end-1) = spike_buffer(:, 2:end);
    spike_buffer(:, end) = 0;
    if ~isempty(fired)
        spike_buffer(fired, end) = 1;
    end

    delayed = find(spike_buffer(:,1));
    if isempty(delayed)
        I_syn = zeros(N,1);
    else
        I_syn = sum(net.W(:, delayed), 2);
    end

    I_tot = I_bg + I_syn;
    dv = 0.04 * v.^2 + 5 * v + 140 - u + I_tot;
    v = v + 0.5 * dv * dt;
    dv = 0.04 * v.^2 + 5 * v + 140 - u + I_tot;
    v = v + 0.5 * dv * dt;
    u = u + net.a .* (net.b .* v - u) * dt;

    v_trace(:, t) = v;
end
end


function data = voltage_to_point_process(v_trace, dt_ms)
N = size(v_trace, 1);
data = struct('times', cell(N,1));

for i = 1:N
    spike_idx = find(v_trace(i, 2:end) >= 30 & v_trace(i, 1:end-1) < 30) + 1;
    data(i).times = (spike_idx(:) * dt_ms) / 1000;
end
end


function [S, f] = cpsd_voltage_spectrum(v_trace, spec)
X = zscore(v_trace, 0, 2).';
[~, f] = cpsd(X(:,1), X(:,1), spec.win, spec.nov, spec.nfft, spec.fs);
nf = numel(f);
n = size(X, 2);
S = zeros(n, n, nf);

for i = 1:n
    for j = i:n
        Pij = cpsd(X(:,i), X(:,j), spec.win, spec.nov, spec.nfft, spec.fs);
        S(i,j,:) = Pij;
        if i ~= j
            S(j,i,:) = conj(Pij);
        end
    end
end
end


function S = stitch_block_spectra(blockSpectra, recset, N)
nf = size(blockSpectra{1}, 3);
S = zeros(N, N, nf);
count = zeros(N, N);

for u = 1:numel(recset)
    idx = recset{u};
    S(idx, idx, :) = S(idx, idx, :) + blockSpectra{u};
    count(idx, idx) = count(idx, idx) + 1;
end

if any(count(:) == 0)
    error('Observation blocks do not cover every spectral entry.');
end

S = S ./ repmat(count, 1, 1, nf);
for k = 1:nf
    S(:,:,k) = (S(:,:,k) + S(:,:,k)') / 2;
end
end


function flow = pdc_to_infoflow(pdc)
[n, ~, nf] = size(pdc);
residual = pdc - repmat(eye(n), 1, 1, nf);
residual_power = min(abs(residual).^2, 1 - 1e-12);
flow = real(-sum(log(1 - residual_power), 3));
% Matches the truth matrix convention: flow(i,j) is j -> i.
flow(1:n+1:end) = 0;
end


function assert_nonempty_spikes(data, label)
spike_counts = cellfun(@numel, {data.times});
if any(spike_counts == 0)
    error(['At least one neuron produced zero spikes during ', label, '. ', ...
           'Increase sim_ms, input_mean, or noise_amp.']);
end
end


function r = offdiag_corr(A, B)
mask = ~eye(size(A, 1));
a = A(mask);
b = B(mask);
good = isfinite(a) & isfinite(b);
a = a(good);
b = b(good);
if numel(a) < 2 || std(a) == 0 || std(b) == 0
    r = NaN;
else
    a = a - mean(a);
    b = b - mean(b);
    r = sum(a .* b) / sqrt(sum(a.^2) * sum(b.^2));
end
end


function plot_flow_matrix(M, labels)
imagesc(M);
axis image;
set(gca, 'XTick', 1:numel(labels), 'XTickLabel', labels, ...
         'YTick', 1:numel(labels), 'YTickLabel', labels, ...
         'TickLabelInterpreter', 'none');
xlabel('Source');
ylabel('Target');
colormap(gca, blue_red_colormap(256));
offdiag = M(~eye(size(M,1)));
cmax = max(offdiag);
if cmax <= 0
    cmax = 1;
end
caxis([0 cmax]);
colorbar;

for i = 1:size(M,1)
    for j = 1:size(M,2)
        if i ~= j
            if M(i,j) > 0.55 * cmax
                txtColor = 'w';
            else
                txtColor = 'k';
            end
            text(j, i, sprintf('%.2f', M(i,j)), ...
                'HorizontalAlignment', 'center', 'Color', txtColor, 'FontSize', 10);
        end
    end
end
end


function cmap = blue_red_colormap(n)
if nargin < 1
    n = 256;
end
blue = [0.08 0.25 0.58];
white = [0.98 0.98 0.96];
red = [0.70 0.08 0.12];
x = linspace(0, 1, n)';
cmap = zeros(n, 3);
for k = 1:3
    cmap(:,k) = interp1([0 0.5 1], [blue(k) white(k) red(k)], x, 'linear');
end
end
