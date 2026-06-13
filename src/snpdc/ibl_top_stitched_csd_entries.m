function T = ibl_top_stitched_csd_entries(result, freqBand, nTop)
%IBL_TOP_STITCHED_CSD_ENTRIES List largest off-diagonal stitched CSD entries.

if nargin < 2 || isempty(freqBand)
    freqIdx = 1:numel(result.freqs);
else
    band = sort(freqBand(:));
    freqIdx = find(result.freqs >= band(1) & result.freqs <= band(end));
end
if nargin < 3 || isempty(nTop)
    nTop = 50;
end

areaNames = string(result.area_names(:));
M = mean(abs(result.S(:,:,freqIdx)), 3, 'omitnan');
M(1:size(M,1)+1:end) = NaN;

[vals, ord] = sort(M(:), 'descend', 'MissingPlacement', 'last');
ord = ord(isfinite(vals));
vals = vals(isfinite(vals));
ord = ord(1:min(nTop, numel(ord)));
vals = vals(1:numel(ord));
[rowIdx, colIdx] = ind2sub(size(M), ord);

T = table(areaNames(rowIdx), areaNames(colIdx), rowIdx, colIdx, vals, ...
    'VariableNames', {'target_area', 'source_area', 'target_index', 'source_index', 'mean_abs_csd'});
end
