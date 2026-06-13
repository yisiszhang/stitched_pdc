function out = ibl_plot_cross_area_csd(cfg, sessionId, bandHz)
%IBL_PLOT_CROSS_AREA_CSD Visualize a saved cross-area CSD summary before stitching.
%
%   out = ibl_plot_cross_area_csd(cfg, sessionId)
%   out = ibl_plot_cross_area_csd(cfg, sessionId, bandHz)
%
% Inputs
%   cfg       output config from ibl_default_config
%   sessionId session identifier matching a file in cfg.cross_spectra_dir
%   bandHz    two-element frequency band, default [1 40]
%
% Output
%   out       struct containing the loaded summary and derived matrices

if nargin < 2 || isempty(sessionId)
    files = dir(fullfile(cfg.cross_spectra_dir, '*.mat'));
    assert(~isempty(files), 'No cross-spectrum summaries found in %s', cfg.cross_spectra_dir);
    [~, ord] = sort({files.name});
    files = files(ord);
    sessionId = erase(files(1).name, '.mat');
end
if nargin < 3 || isempty(bandHz)
    bandHz = [1 40];
end

file = fullfile(cfg.cross_spectra_dir, [char(sessionId) '.mat']);
assert(exist(file, 'file') == 2, 'Cross-spectrum summary not found: %s', file);

tmp = load(file, 'summary');
summary = tmp.summary;
freqs = summary.freqs(:);
mask = freqs >= bandHz(1) & freqs <= min(bandHz(2), freqs(end));
assert(any(mask), 'No frequencies fall inside [%g %g] Hz', bandHz(1), bandHz(2));

S = summary.cross_spectrum;
auto = summary.auto_spectrum;
areaNames = string(summary.area_names(:));
nAreas = numel(areaNames);

bandAbsCsd = mean(abs(S(:,:,mask)), 3, 'omitnan');
bandRealCsd = mean(real(S(:,:,mask)), 3, 'omitnan');
bandImagCsd = mean(imag(S(:,:,mask)), 3, 'omitnan');

coh = nan(nAreas, nAreas, sum(mask));
fIdx = find(mask);
for k = 1:numel(fIdx)
    idx = fIdx(k);
    sdiag = real(diag(S(:,:,idx)));
    denom = sqrt(max(sdiag, 0) * max(sdiag, 0).');
    denom(denom == 0) = nan;
    coh(:,:,k) = abs(S(:,:,idx)) ./ denom;
end
bandCoh = mean(coh, 3, 'omitnan');
bandCoh(1:nAreas+1:end) = 1;

figure('Name', sprintf('IBL CSD: %s', summary.session_id), 'Color', 'w');
tiledlayout(2,2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
imagesc(bandAbsCsd);
axis image;
colorbar;
title(sprintf('|CSD| mean, %g-%g Hz', bandHz(1), bandHz(2)));
set(gca, 'XTick', 1:nAreas, 'XTickLabel', areaNames, ...
    'YTick', 1:nAreas, 'YTickLabel', areaNames, 'XTickLabelRotation', 90);

nexttile;
imagesc(bandCoh, [0 1]);
axis image;
colorbar;
title(sprintf('Coherence mean, %g-%g Hz', bandHz(1), bandHz(2)));
set(gca, 'XTick', 1:nAreas, 'XTickLabel', areaNames, ...
    'YTick', 1:nAreas, 'YTickLabel', areaNames, 'XTickLabelRotation', 90);

nexttile;
imagesc(bandRealCsd);
axis image;
colorbar;
title(sprintf('Real(CSD) mean, %g-%g Hz', bandHz(1), bandHz(2)));
set(gca, 'XTick', 1:nAreas, 'XTickLabel', areaNames, ...
    'YTick', 1:nAreas, 'YTickLabel', areaNames, 'XTickLabelRotation', 90);

nexttile;
imagesc(bandImagCsd);
axis image;
colorbar;
title(sprintf('Imag(CSD) mean, %g-%g Hz', bandHz(1), bandHz(2)));
set(gca, 'XTick', 1:nAreas, 'XTickLabel', areaNames, ...
    'YTick', 1:nAreas, 'YTickLabel', areaNames, 'XTickLabelRotation', 90);

sgtitle(sprintf('Cross-Area CSD Summary: %s', summary.session_id), 'Interpreter', 'none');

out.summary = summary;
out.band = bandHz;
out.band_mask = mask;
out.band_abs_csd = bandAbsCsd;
out.band_real_csd = bandRealCsd;
out.band_imag_csd = bandImagCsd;
out.band_coherence = bandCoh;
out.file = file;
end
