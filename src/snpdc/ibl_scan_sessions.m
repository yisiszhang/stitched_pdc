function scan = ibl_scan_sessions(cfg)
%IBL_SCAN_SESSIONS Discover and filter IBL spontaneous sessions from local cache.

regionTable = ibl_load_region_table(cfg.cache_root);
files = dir(fullfile(cfg.cache_root, '*', 'Subjects', '*', '*', '*', 'alf', '_ibl_passivePeriods.intervalsTable.csv'));
files = files(~startsWith({files.name}, '._'));

if isfield(cfg, 'verbose') && cfg.verbose
    fprintf('[ibl_scan_sessions] Found %d passive-period files under %s\n', numel(files), cfg.cache_root);
end

sessionRows = struct('session_id', {}, 'session_root', {}, 'lab', {}, 'subject', {}, ...
    'date', {}, 'number', {}, 'sp_start', {}, 'sp_end', {}, 'sp_dur', {}, 'areas', {}, 'n_areas', {});
areaCount = containers.Map('KeyType', 'char', 'ValueType', 'double');
scanTimer = tic;

for k = 1:numel(files)
    csvFile = fullfile(files(k).folder, files(k).name);
    sessionRoot = fileparts(files(k).folder);
    parts = strsplit(sessionRoot, filesep);
    number = parts{end};
    date = parts{end-1};
    subject = parts{end-2};
    lab = parts{end-4};

    interval = ibl_read_spontaneous_interval(csvFile);
    if isempty(interval)
        continue;
    end

    spDur = interval(2) - interval(1);
    if spDur < cfg.min_sp_dur
        local_progress(cfg, k, numel(files), sessionRoot, scanTimer, 0, 'short_spontaneous');
        continue;
    end

    try
        areaStruct = ibl_load_area_spikes(sessionRoot, regionTable, cfg);
    catch ME
        warning('Skipping %s: %s', sessionRoot, ME.message);
        local_progress(cfg, k, numel(files), sessionRoot, scanTimer, 0, 'error');
        continue;
    end

    if numel(areaStruct) < 2
        local_progress(cfg, k, numel(files), sessionRoot, scanTimer, numel(areaStruct), 'too_few_areas');
        continue;
    end

    areas = sort(string({areaStruct.area}));
    for j = 1:numel(areas)
        key = char(areas(j));
        if ~areaCount.isKey(key)
            areaCount(key) = 1;
        else
            areaCount(key) = areaCount(key) + 1;
        end
    end

    row.session_id = sprintf('%s__%s__%s__%s', lab, subject, date, number);
    row.session_root = sessionRoot;
    row.lab = lab;
    row.subject = subject;
    row.date = date;
    row.number = number;
    row.sp_start = interval(1);
    row.sp_end = interval(2);
    row.sp_dur = spDur;
    row.areas = areas;
    row.n_areas = numel(areas);
    sessionRows(end+1) = row; %#ok<AGROW>
    local_progress(cfg, k, numel(files), sessionRoot, scanTimer, numel(areas), 'kept');
end

allAreas = sort(string(keys(areaCount)));
if isempty(allAreas)
    error('No qualifying sessions were found under %s', cfg.cache_root);
end

coObs = false(numel(allAreas));
for k = 1:numel(sessionRows)
    [~, idx] = ismember(sessionRows(k).areas, allAreas);
    idx = idx(idx > 0);
    coObs(idx, idx) = true;
end

components = local_connected_components(coObs, allAreas);
largest = components{1};

qualifying = sessionRows(arrayfun(@(s) sum(ismember(s.areas, largest)) >= 2, sessionRows));

scan.sessions = sessionRows;
scan.qualifying_sessions = qualifying;
scan.area_names = allAreas;
scan.co_observation = coObs;
scan.components = components;
scan.largest_component = largest;
scan.region_table = regionTable;

save(cfg.scan_file, 'scan', '-v7.3');
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


function local_progress(cfg, idx, total, sessionRoot, scanTimer, nAreas, status)
if ~isfield(cfg, 'verbose') || ~cfg.verbose
    return;
end

every = 10;
if isfield(cfg, 'progress_every') && ~isempty(cfg.progress_every)
    every = cfg.progress_every;
end

if ~(idx == 1 || idx == total || strcmp(status, 'kept') || mod(idx, every) == 0)
    return;
end

elapsed = toc(scanTimer);
rate = elapsed / max(idx, 1);
remaining = rate * max(total - idx, 0);

[~, sessionName] = fileparts(sessionRoot);
fprintf('[ibl_scan_sessions] %4d/%4d  %-14s  areas=%2d  elapsed=%6.1fs  eta=%6.1fs  %s\n', ...
    idx, total, status, nAreas, elapsed, remaining, sessionName);
end
