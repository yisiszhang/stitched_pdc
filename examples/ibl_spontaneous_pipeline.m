% ibl_spontaneous_pipeline.m
% Local IBL spontaneous Neuropixels pipeline:
%   1. scan/filter sessions from ONE cache
%   2. compute within-area spectra/PCA using Chronux point-process multitaper
%   3. apply PCA QC and recompute the largest connected component
%   4. compute cross-area CSD only within the passed component
%   5. stitch area-level spectra across sessions/animals
%   6. estimate stitched nonparametric PDC

clear; close all; clc;
setup;

cfg = ibl_default_config( ...
    '/Volumes/Extreme SSD/data/neuropixel/ONE/openalyx.internationalbrainlab.org', ...
    fullfile(pwd, 'ibl_output'));

scan = ibl_scan_sessions(cfg);
areaPca = ibl_compute_area_pca(scan, cfg); %#ok<NASGU>
qc = ibl_build_pca_qc_graph(cfg, scan);
crossSummaries = ibl_compute_cross_area_csd(qc, cfg); %#ok<NASGU>
result = ibl_stitch_saved_spectra(cfg, qc);

disp('Largest connected component areas:');
disp(qc.largest_component');

figure('Color', 'w', 'Position', [80 80 900 700]);
imagesc(qc.co_observation);
axis image;
set(gca, 'XTick', 1:numel(qc.area_names), 'XTickLabel', qc.area_names, ...
    'YTick', 1:numel(qc.area_names), 'YTickLabel', qc.area_names, ...
    'XTickLabelRotation', 60, 'TickLabelInterpreter', 'none');
title('Area co-observation graph after PCA QC');
colorbar;

fMask = result.freqs >= 1 & result.freqs <= min(80, result.freqs(end));
meanPdc = mean(abs(result.PDC(:,:,fMask)), 3);
meanPdc(1:size(meanPdc,1)+1:end) = 0;

figure('Color', 'w', 'Position', [120 120 900 700]);
imagesc(meanPdc);
axis image;
set(gca, 'XTick', 1:numel(result.area_names), 'XTickLabel', result.area_names, ...
    'YTick', 1:numel(result.area_names), 'YTickLabel', result.area_names, ...
    'XTickLabelRotation', 60, 'TickLabelInterpreter', 'none');
xlabel('Source area');
ylabel('Target area');
title('Mean stitched PDC magnitude (1-80 Hz)');
colorbar;

infoflow = ibl_pdc_to_infoflow(result.PDC(:,:,fMask));
infoflow(1:size(infoflow,1)+1:end) = 0;

figure('Color', 'w', 'Position', [160 160 900 700]);
imagesc(infoflow);
axis image;
set(gca, 'XTick', 1:numel(result.area_names), 'XTickLabel', result.area_names, ...
    'YTick', 1:numel(result.area_names), 'YTickLabel', result.area_names, ...
    'XTickLabelRotation', 60, 'TickLabelInterpreter', 'none');
xlabel('Source area');
ylabel('Target area');
title('Integrated information flow from nonparametric PDC (1-80 Hz)');
colorbar;
