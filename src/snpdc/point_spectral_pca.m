function summary = point_spectral_pca(spikeTimes, duration, params)
%POINT_SPECTRAL_PCA Spectral PCA for a population of point processes.
%
%   summary = point_spectral_pca(spikeTimes, duration, params)
%
% Inputs
%   spikeTimes : cell array; spikeTimes{i} is the spike-time vector of unit i
%   duration   : recording duration in seconds
%   params     : Chronux-style struct. Required fields are those accepted by
%                CrossSpecMatpt. Optional fields:
%                   n_pcs          default 3
%                   pc_band        default [1 100]
%                   target_n_freqs default [] (no reduction)
%
% Output
%   summary.Saa                  neuron-by-neuron point-process spectrum
%   summary.freqs                frequency grid
%   summary.eigenvalues          spectral PCA eigenvalues
%   summary.eigenvectors         spectral PCA loadings
%   summary.explained_variance   frequency-wise explained variance
%   summary.fixed_loading        band-averaged PC1 loading
%   summary.fixed_loading_eigval band-averaged PC1 eigenvalue

assert(exist('CrossSpecMatpt', 'file') == 2, ...
    'Chronux CrossSpecMatpt not found. Add Chronux to the MATLAB path.');

if ~isfield(params, 'n_pcs'); params.n_pcs = 3; end
if ~isfield(params, 'pc_band'); params.pc_band = [1 100]; end
if ~isfield(params, 'target_n_freqs'); params.target_n_freqs = []; end

data = local_times_to_struct(spikeTimes);
[pxy, ~, ~, ~, ~, f] = CrossSpecMatpt(data, params.win, duration, params);
Saa = permute(pxy, [3, 2, 1]);
for i = 1:size(Saa,3)
    Saa(:,:,i) = (Saa(:,:,i) + Saa(:,:,i)') / 2;
end
Saa(~isfinite(Saa)) = 0;

if ~isempty(params.target_n_freqs)
    [Saa, f] = reduce_frequency_grid(Saa, f, params.target_n_freqs);
end

[evals, evecs, explVar] = spectral_pca(Saa, params.n_pcs);
[vfix, dfix, Sbar] = spectral_pc_loading(Saa, f, params.pc_band);

summary = struct();
summary.Saa = Saa;
summary.freqs = f(:);
summary.eigenvalues = evals;
summary.eigenvectors = evecs;
summary.explained_variance = explVar;
summary.fixed_loading = vfix;
summary.fixed_loading_eigval = dfix;
summary.fixed_loading_band_spectrum = Sbar;
summary.pc_band = params.pc_band;
summary.sign_convention = "largest_abs_loading_real_positive";
end


function data = local_times_to_struct(timesCell)
n = numel(timesCell);
data = struct('times', cell(n,1));
for i = 1:n
    data(i).times = timesCell{i}(:);
end
end
