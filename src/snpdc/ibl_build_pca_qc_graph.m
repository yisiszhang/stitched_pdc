function qc = ibl_build_pca_qc_graph(cfg, scan)
%IBL_BUILD_PCA_QC_GRAPH Build area co-observation graph after area-level PCA QC.

files = dir(fullfile(cfg.area_pca_dir, '*.mat'));
assert(~isempty(files), 'No area PCA files found in %s', cfg.area_pca_dir);

sessionRows = struct('session_id', {}, 'area_names', {}, 'file', {}, 'n_areas', {});
areaCount = containers.Map('KeyType', 'char', 'ValueType', 'double');

for k = 1:numel(files)
    tmp = load(fullfile(files(k).folder, files(k).name), 'summary');
    summary = tmp.summary;
    keep = summary.pass_qc_area(:);
    areaNames = string(summary.area_names(keep));
    if numel(areaNames) < 2
        continue;
    end
    areaNames = sort(areaNames);

    row.session_id = summary.session_id;
    row.area_names = areaNames;
    row.file = fullfile(files(k).folder, files(k).name);
    row.n_areas = numel(areaNames);
    sessionRows(end+1) = row; %#ok<AGROW>

    for j = 1:numel(areaNames)
        key = char(areaNames(j));
        if ~areaCount.isKey(key)
            areaCount(key) = 1;
        else
            areaCount(key) = areaCount(key) + 1;
        end
    end
end

allAreas = sort(string(keys(areaCount)));
if isempty(allAreas)
    error('No areas passed PCA QC in at least one session.');
end

minSupport = 1;
if isfield(cfg, 'min_sessions_per_area') && ~isempty(cfg.min_sessions_per_area)
    minSupport = cfg.min_sessions_per_area;
end

supportCounts = zeros(numel(allAreas),1);
for i = 1:numel(allAreas)
    supportCounts(i) = areaCount(char(allAreas(i)));
end

supportMask = supportCounts >= minSupport;
filteredAreas = allAreas(supportMask);
if isempty(filteredAreas)
    error('No areas met min_sessions_per_area=%d after PCA QC.', minSupport);
end

if cfg.verbose
    fprintf('[ibl_build_pca_qc_graph] Areas passing PCA QC: %d\n', numel(allAreas));
    fprintf('[ibl_build_pca_qc_graph] Areas kept with min_sessions_per_area=%d: %d\n', ...
        minSupport, numel(filteredAreas));
end

coObs = false(numel(allAreas));
for k = 1:numel(sessionRows)
    [~, idx] = ismember(sessionRows(k).area_names, allAreas);
    idx = idx(idx > 0);
    coObs(idx, idx) = true;
end

coObs = coObs(supportMask, supportMask);
allAreas = filteredAreas;

components = local_connected_components(coObs, allAreas);
largest = strings(0,1);
if ~isempty(components)
    largest = components{1};
end

qualifying = sessionRows(arrayfun(@(s) sum(ismember(s.area_names, largest)) >= 2, sessionRows));

qc.sessions = sessionRows;
qc.qualifying_sessions = qualifying;
qc.area_names = allAreas;
qc.co_observation = coObs;
qc.components = components;
qc.largest_component = largest;
qc.area_support_count = supportCounts(supportMask);
qc.min_sessions_per_area = minSupport;
qc.scan = scan;

save(cfg.pca_qc_file, 'qc', '-v7.3');
end


function components = local_connected_components(adj, labels)
n = size(adj,1);
seen = false(n,1);
components = {};
for i = 1:n
    if seen(i)
        continue;
    end
    queue = i;
    comp = [];
    seen(i) = true;
    while ~isempty(queue)
        cur = queue(1);
        queue(1) = [];
        comp(end+1) = cur; %#ok<AGROW>
        nb = find(adj(cur,:));
        nb = nb(~seen(nb));
        seen(nb) = true;
        queue = [queue nb]; %#ok<AGROW>
    end
    components{end+1} = sort(labels(comp)); %#ok<AGROW>
end
[~, ord] = sort(cellfun(@numel, components), 'descend');
components = components(ord);
end
