function fp = ibl_compare_cross_session_spectral_fingerprints(crossCand, cfg, varargin)
%IBL_COMPARE_CROSS_SESSION_SPECTRAL_FINGERPRINTS Compare area latent spectra.
%
%   fp = ibl_compare_cross_session_spectral_fingerprints(crossCand, cfg)
%
% This diagnostic compares source and target sessions before CSD/PDC. It
% uses cached area-PCA summaries and dimension-invariant spectral features:
% PC explained variance by band, effective rank, eigenvalue spectra, firing
% rate summaries, and neuron count.

p = inputParser;
p.addParameter('CandidateIndex', 1, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('Bands', [1 10; 10 30; 30 80], @(x)isnumeric(x)&&size(x,2)==2);
p.addParameter('BandNames', strings(0,1), @(x)isstring(x)||iscellstr(x)||isempty(x));
p.addParameter('MaxPCs', 3, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('MakeFigure', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('OutputFile', "", @(x)isstring(x)||ischar(x));
p.parse(varargin{:});
opt = p.Results;
opt.BandNames = string(opt.BandNames);
if isempty(opt.BandNames)
    opt.BandNames = arrayfun(@(i) sprintf('%g-%gHz', opt.Bands(i,1), opt.Bands(i,2)), ...
        (1:size(opt.Bands,1)).', 'UniformOutput', false);
    opt.BandNames = string(opt.BandNames);
end

if isfield(crossCand, 'sets')
    cand = crossCand.sets(opt.CandidateIndex);
else
    cand = crossCand(opt.CandidateIndex);
end

targetSummary = local_load_summary(cand.target.pca_file);
targetAreas = cand.covered_area_names(:);
targetFeatures = local_session_features(targetSummary, targetAreas, opt);

rows = struct('area', {}, 'source_session_id', {}, 'source_lab', {}, ...
    'source_subject', {}, 'target_session_id', {}, ...
    'feature_corr', {}, 'explained_var_corr', {}, 'eig_spectrum_corr', {}, ...
    'power_corr', {}, 'effective_rank_diff', {}, 'pc1_band_diff', {}, ...
    'n_neurons_target', {}, 'n_neurons_source', {}, ...
    'mean_fr_target', {}, 'mean_fr_source', {});

for s = 1:numel(cand.sources)
    src = cand.sources(s);
    srcSummary = local_load_summary(src.pca_file);
    sharedAreas = intersect(string(src.area_names(:)), targetAreas, 'stable');
    srcFeatures = local_session_features(srcSummary, sharedAreas, opt);
    targetShared = local_session_features(targetSummary, sharedAreas, opt);
    for a = 1:numel(sharedAreas)
        row.area = sharedAreas(a);
        row.source_session_id = src.session_id;
        row.source_lab = src.lab;
        row.source_subject = src.subject;
        row.target_session_id = cand.target.session_id;
        row.feature_corr = local_corr(srcFeatures(a).feature_vector, targetShared(a).feature_vector);
        row.explained_var_corr = local_corr(srcFeatures(a).explained_vector, targetShared(a).explained_vector);
        row.eig_spectrum_corr = local_corr(srcFeatures(a).eig_vector, targetShared(a).eig_vector);
        row.power_corr = local_corr(srcFeatures(a).power_vector, targetShared(a).power_vector);
        row.effective_rank_diff = srcFeatures(a).effective_rank_mean - targetShared(a).effective_rank_mean;
        row.pc1_band_diff = srcFeatures(a).pc1_band_mean - targetShared(a).pc1_band_mean;
        row.n_neurons_target = targetShared(a).n_neurons;
        row.n_neurons_source = srcFeatures(a).n_neurons;
        row.mean_fr_target = targetShared(a).mean_fr;
        row.mean_fr_source = srcFeatures(a).mean_fr;
        rows(end+1) = row; %#ok<AGROW>
    end
end

if isempty(rows)
    areaTable = table();
else
    areaTable = struct2table(rows);
    areaTable = sortrows(areaTable, {'area', 'feature_corr'}, {'ascend', 'descend'});
end

areaSummary = local_area_summary(areaTable);

fp.candidate_index = opt.CandidateIndex;
fp.target_session_id = cand.target.session_id;
fp.target_area_names = targetAreas;
fp.target_features = targetFeatures;
fp.area_table = areaTable;
fp.area_summary = areaSummary;
fp.options = opt;
fp.created_at = string(datetime('now'));

if strlength(string(opt.OutputFile)) > 0
    save(opt.OutputFile, 'fp', '-v7.3');
end

if opt.MakeFigure
    local_plot_fingerprints(fp);
end
end


function summary = local_load_summary(pcaFile)
tmp = load(pcaFile, 'summary');
summary = tmp.summary;
end


function features = local_session_features(summary, areaNames, opt)
summaryAreas = string(summary.area_names(:));
features = struct('area', {}, 'feature_vector', {}, 'explained_vector', {}, ...
    'eig_vector', {}, 'power_vector', {}, 'effective_rank_mean', {}, ...
    'pc1_band_mean', {}, 'n_neurons', {}, 'mean_fr', {});

for a = 1:numel(areaNames)
    area = areaNames(a);
    idx = find(summaryAreas == area, 1);
    assert(~isempty(idx), 'Area %s not found in %s.', area, summary.session_id);

    freqs = summary.freqs(:);
    expl = summary.expl_var{idx};
    eigvals = summary.eigenvalues{idx};
    maxPcs = min([opt.MaxPCs, size(expl,1), size(eigvals,1)]);

    explainedByBand = local_band_means(expl(1:maxPcs,:), freqs, opt.Bands);
    eigByBand = local_band_means(log1p(max(real(eigvals(1:maxPcs,:)), 0)), freqs, opt.Bands);
    powerByBand = local_band_means(log1p(max(real(summary.auto_spectrum(idx,:)), 0)), freqs, opt.Bands);
    effectiveRank = local_effective_rank(eigvals, freqs, opt.Bands);

    fr = summary.firing_rates{idx};
    if isempty(fr)
        meanFr = NaN;
    else
        meanFr = mean(fr, 'omitnan');
    end

    featureVector = [explainedByBand(:); eigByBand(:); powerByBand(:); effectiveRank(:); ...
        log1p(summary.n_neurons(idx)); log1p(meanFr)];

    feat.area = area;
    feat.feature_vector = featureVector;
    feat.explained_vector = explainedByBand(:);
    feat.eig_vector = eigByBand(:);
    feat.power_vector = powerByBand(:);
    feat.effective_rank_mean = mean(effectiveRank, 'omitnan');
    feat.pc1_band_mean = mean(explainedByBand(1,:), 'omitnan');
    feat.n_neurons = summary.n_neurons(idx);
    feat.mean_fr = meanFr;
    features(end+1) = feat; %#ok<AGROW>
end
end


function out = local_band_means(values, freqs, bands)
if isrow(values)
    values = values(:).';
end
nRows = size(values, 1);
nBands = size(bands, 1);
out = nan(nRows, nBands);
for b = 1:nBands
    mask = freqs >= bands(b,1) & freqs < bands(b,2);
    if ~any(mask)
        mask = freqs >= bands(b,1) & freqs <= bands(b,2);
    end
    if any(mask)
        out(:,b) = mean(values(:,mask), 2, 'omitnan');
    end
end
end


function er = local_effective_rank(eigvals, freqs, bands)
eigvals = max(real(eigvals), 0);
nBands = size(bands, 1);
er = nan(1, nBands);
for b = 1:nBands
    mask = freqs >= bands(b,1) & freqs < bands(b,2);
    if ~any(mask)
        mask = freqs >= bands(b,1) & freqs <= bands(b,2);
    end
    if ~any(mask)
        continue;
    end
    vals = eigvals(:,mask);
    vals = vals ./ max(sum(vals, 1), eps);
    entropy = -sum(vals .* log(max(vals, eps)), 1);
    er(b) = mean(exp(entropy), 'omitnan');
end
end


function r = local_corr(a, b)
a = real(a(:));
b = real(b(:));
good = isfinite(a) & isfinite(b);
if nnz(good) < 3 || std(a(good)) == 0 || std(b(good)) == 0
    r = NaN;
else
    r = corr(a(good), b(good), 'type', 'Pearson');
end
end


function T = local_area_summary(areaTable)
if isempty(areaTable)
    T = table();
    return;
end
areas = unique(areaTable.area, 'stable');
rows = struct('area', {}, 'n_sources', {}, ...
    'feature_corr_median', {}, 'feature_corr_max', {}, ...
    'explained_var_corr_median', {}, 'eig_spectrum_corr_median', {}, ...
    'power_corr_median', {}, 'effective_rank_absdiff_median', {}, ...
    'pc1_band_absdiff_median', {});
for i = 1:numel(areas)
    mask = areaTable.area == areas(i);
    row.area = areas(i);
    row.n_sources = nnz(mask);
    row.feature_corr_median = local_nanmedian(areaTable.feature_corr(mask));
    row.feature_corr_max = local_nanmax(areaTable.feature_corr(mask));
    row.explained_var_corr_median = local_nanmedian(areaTable.explained_var_corr(mask));
    row.eig_spectrum_corr_median = local_nanmedian(areaTable.eig_spectrum_corr(mask));
    row.power_corr_median = local_nanmedian(areaTable.power_corr(mask));
    row.effective_rank_absdiff_median = local_nanmedian(abs(areaTable.effective_rank_diff(mask)));
    row.pc1_band_absdiff_median = local_nanmedian(abs(areaTable.pc1_band_diff(mask)));
    rows(end+1) = row; %#ok<AGROW>
end
T = struct2table(rows);
T = sortrows(T, 'feature_corr_median', 'descend');
end


function y = local_nanmedian(x)
x = x(isfinite(x));
if isempty(x)
    y = NaN;
else
    y = median(x);
end
end


function y = local_nanmax(x)
x = x(isfinite(x));
if isempty(x)
    y = NaN;
else
    y = max(x);
end
end


function local_plot_fingerprints(fp)
T = fp.area_table;
if isempty(T)
    return;
end

figure('Color', 'w', 'Position', [100 100 1180 430]);
tiledlayout(1,3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
boxchart(categorical(T.area), T.feature_corr);
ylabel('Feature correlation');
title('Source-target spectral fingerprint similarity');
grid on;
set(gca, 'XTickLabelRotation', 45, 'TickLabelInterpreter', 'none');

nexttile;
scatter(T.explained_var_corr, T.eig_spectrum_corr, 45, 'filled', 'MarkerFaceAlpha', 0.7);
xlabel('Explained-variance corr');
ylabel('Eigen-spectrum corr');
title('Feature components');
grid on;

nexttile;
scatter(T.feature_corr, abs(T.pc1_band_diff), 45, 'filled', 'MarkerFaceAlpha', 0.7);
xlabel('Feature correlation');
ylabel('|PC1 explained diff|');
title('PC1 similarity vs full fingerprint');
grid on;
end
