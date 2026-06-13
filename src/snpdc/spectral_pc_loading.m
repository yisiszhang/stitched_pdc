function [vfix, dfix, Sbar] = spectral_pc_loading(Saa, freqs, band)
%SPECTRAL_PC_LOADING Fixed loading from a band-averaged spectrum.
%
%   [vfix, dfix, Sbar] = spectral_pc_loading(Saa, freqs, band)
%
% This is useful when projecting point-process activity from multiple
% neurons in one area to a single area-level latent signal before estimating
% cross-area spectra.

if nargin < 3 || isempty(band)
    band = [min(freqs(:)) max(freqs(:))];
end

mask = freqs >= band(1) & freqs <= min(band(2), freqs(end));
if ~any(mask)
    mask = true(size(freqs));
end

Sbar = mean(Saa(:,:,mask), 3, 'omitnan');
Sbar(~isfinite(Sbar)) = 0;
Sbar = (Sbar + Sbar') / 2;

[V, D] = eig(Sbar, 'vector');
[d, ord] = sort(real(D), 'descend');
d = max(d, 0);
V = V(:, ord);

vfix = orient_loading(V(:,1));
dfix = d(1);

if ~any(isfinite(vfix)) || norm(vfix) == 0
    vfix = zeros(size(vfix));
    vfix(1) = 1;
else
    vfix = vfix / norm(vfix);
    vfix = orient_loading(vfix);
end
end
