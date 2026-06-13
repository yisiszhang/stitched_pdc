function out = ibl_inspect_stitched_csd_entry(result, cfg, areaA, areaB, freqOrBand)
%IBL_INSPECT_STITCHED_CSD_ENTRY Trace one stitched CSD entry to source blocks.
%
%   out = ibl_inspect_stitched_csd_entry(result, cfg, areaA, areaB, freqOrBand)
%
% areaA/areaB can be acronyms or numeric indices in result.area_names.
% freqOrBand can be a scalar frequency/index or a two-element Hz band.

areaNames = string(result.area_names(:));
idxA = local_area_index(areaA, areaNames);
idxB = local_area_index(areaB, areaNames);
freqIdx = local_frequency_indices(freqOrBand, result.freqs(:));

files = dir(fullfile(cfg.cross_spectra_dir, '*.mat'));
assert(~isempty(files), 'No cross-spectra files found in %s.', cfg.cross_spectra_dir);

vals = [];
rows = struct('session_id', {}, 'value', {}, 'magnitude', {}, ...
    'auto_a', {}, 'auto_b', {}, 'coherence', {}, 'n_freqs', {});

for k = 1:numel(files)
    tmp = load(fullfile(files(k).folder, files(k).name), 'summary');
    s = tmp.summary;
    names = string(s.area_names(:));
    ia = find(names == areaNames(idxA), 1);
    ib = find(names == areaNames(idxB), 1);
    if isempty(ia) || isempty(ib)
        continue;
    end
    assert(numel(s.freqs) == numel(result.freqs) && max(abs(s.freqs(:) - result.freqs(:))) < 1e-10, ...
        'Frequency mismatch in %s.', s.session_id);

    v = squeeze(mean(s.cross_spectrum(ia, ib, freqIdx), 3, 'omitnan'));
    aa = squeeze(mean(s.cross_spectrum(ia, ia, freqIdx), 3, 'omitnan'));
    bb = squeeze(mean(s.cross_spectrum(ib, ib, freqIdx), 3, 'omitnan'));
    coh = abs(v) ./ sqrt(max(real(aa),0) .* max(real(bb),0));

    row.session_id = string(s.session_id);
    row.value = v;
    row.magnitude = abs(v);
    row.auto_a = real(aa);
    row.auto_b = real(bb);
    row.coherence = coh;
    row.n_freqs = numel(freqIdx);
    rows(end+1) = row; %#ok<AGROW>
    vals(end+1,1) = v; %#ok<AGROW>
end

completedVal = squeeze(mean(result.S(idxA, idxB, freqIdx), 3, 'omitnan'));
if isfield(result, 'Sreg')
    regularizedVal = squeeze(mean(result.Sreg(idxA, idxB, freqIdx), 3, 'omitnan'));
else
    regularizedVal = NaN;
end

if isempty(vals)
    observedMean = NaN;
    observedMedian = NaN;
    observedN = 0;
else
    observedMean = mean(vals, 'omitnan');
    observedMedian = median(vals, 'omitnan');
    observedN = numel(vals);
end

T = struct2table(rows);
if ~isempty(T)
    T = sortrows(T, 'magnitude', 'descend');
end

out.area_a = areaNames(idxA);
out.area_b = areaNames(idxB);
out.index_a = idxA;
out.index_b = idxB;
out.freqs = result.freqs(freqIdx);
out.frequency_indices = freqIdx;
out.n_observed_sessions = observedN;
out.observed_mean = observedMean;
out.observed_median = observedMedian;
out.completed_value = completedVal;
out.regularized_value = regularizedVal;
out.completion_delta_from_observed_mean = completedVal - observedMean;
out.completion_delta_magnitude = abs(completedVal - observedMean);
out.session_table = T;

fprintf('[ibl_inspect_stitched_csd_entry] %s -> %s, %d freq bins\n', ...
    out.area_b, out.area_a, numel(freqIdx));
fprintf('  observed sessions: %d\n', observedN);
fprintf('  observed mean magnitude: %.4g\n', abs(observedMean));
fprintf('  completed magnitude:     %.4g\n', abs(completedVal));
fprintf('  |completed - observed mean|: %.4g\n', out.completion_delta_magnitude);
if observedN == 0
    fprintf('  interpretation: pair was not directly observed; value is from completion.\n');
else
    fprintf('  interpretation: pair was directly observed; compare session_table to completed value.\n');
end
end


function idx = local_area_index(area, areaNames)
if isnumeric(area)
    idx = area;
else
    idx = find(areaNames == string(area), 1);
end
assert(~isempty(idx) && idx >= 1 && idx <= numel(areaNames), 'Area not found: %s', string(area));
end


function idx = local_frequency_indices(freqOrBand, freqs)
if nargin < 1 || isempty(freqOrBand)
    idx = 1:numel(freqs);
    return;
end
if isscalar(freqOrBand)
    if freqOrBand >= 1 && freqOrBand <= numel(freqs) && abs(freqOrBand - round(freqOrBand)) < eps
        idx = round(freqOrBand);
    else
        [~, idx] = min(abs(freqs - freqOrBand));
    end
else
    band = sort(freqOrBand(:));
    idx = find(freqs >= band(1) & freqs <= band(end));
end
assert(~isempty(idx), 'No frequencies selected.');
idx = idx(:)';
end
