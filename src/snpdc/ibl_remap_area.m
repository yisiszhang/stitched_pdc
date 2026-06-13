function parent = ibl_remap_area(acronym, areaMap)
%IBL_REMAP_AREA Map granular Allen labels to parent areas.

if nargin < 2 || isempty(areaMap)
    areaMap = ibl_area_map();
end

if areaMap.map.isKey(acronym)
    parent = areaMap.map(acronym);
else
    parent = '';
end

if any(strcmp(parent, areaMap.exclude))
    parent = '';
    return;
end

if strlength(string(parent)) == 0
    parent = '';
    return;
end

if ~any(areaMap.allowed_outputs == string(parent))
    parent = '';
end
end
