function outFile = ibl_write_area_list(areaNames, outFile)
%IBL_WRITE_AREA_LIST Write one acronym per line for downstream Allen extraction.

if nargin < 2 || isempty(outFile)
    outFile = fullfile(pwd, 'ibl_output', 'allen_area_list.txt');
end

areaNames = string(areaNames(:));
fid = fopen(outFile, 'w');
assert(fid > 0, 'Could not open %s for writing.', outFile);
cleanup = onCleanup(@() fclose(fid));
for i = 1:numel(areaNames)
    fprintf(fid, '%s\n', areaNames(i));
end
end
