function areaStruct = ibl_load_area_spikes(sessionRoot, regionTable, cfg)
%IBL_LOAD_AREA_SPIKES Load spontaneous spike times grouped by remapped area.

interval = ibl_read_spontaneous_interval(fullfile(sessionRoot, 'alf', '_ibl_passivePeriods.intervalsTable.csv'));
assert(~isempty(interval), 'Could not read spontaneous interval for %s', sessionRoot);
spStart = interval(1);
spEnd = interval(2);
spDur = spEnd - spStart;

probeDirs = dir(fullfile(sessionRoot, 'alf', 'probe*', 'pykilosort'));
areaBuckets = containers.Map('KeyType', 'char', 'ValueType', 'any');

for p = 1:numel(probeDirs)
    baseDir = fullfile(probeDirs(p).folder, probeDirs(p).name);
    spikesTimesFile = fullfile(baseDir, 'spikes.times.npy');
    spikesClusFile = fullfile(baseDir, 'spikes.clusters.npy');
    clusterChanFile = fullfile(baseDir, 'clusters.channels.npy');
    brainIdsFile = fullfile(baseDir, 'channels.brainLocationIds_ccf_2017.npy');

    if ~all(cellfun(@(f) exist(f, 'file') == 2, ...
            {spikesTimesFile, spikesClusFile, clusterChanFile, brainIdsFile}))
        continue;
    end

    metricsFile = fullfile(baseDir, 'clusters.metrics.pqt');
    goodMask = ibl_load_good_clusters(metricsFile);

    spikeTimes = double(ibl_read_npy(spikesTimesFile));
    spikeClus = double(ibl_read_npy(spikesClusFile));
    cluChan = double(ibl_read_npy(clusterChanFile));
    brainIds = double(ibl_read_npy(brainIdsFile));

    nClusters = numel(cluChan);
    if isempty(goodMask)
        goodMask = true(nClusters, 1);
    else
        goodMask = logical(goodMask(:));
        if numel(goodMask) > nClusters
            goodMask = goodMask(1:nClusters);
        end
        if numel(goodMask) < nClusters
            goodMask(end+1:nClusters) = false;
        end
    end

    clusterIds = find(goodMask) - 1;  % IBL cluster ids are 0-based in spikes.clusters
    for c = 1:numel(clusterIds)
        cid = clusterIds(c);
        chanIdx = cluChan(cid + 1) + 1;
        if chanIdx < 1 || chanIdx > numel(brainIds)
            continue;
        end

        regionId = brainIds(chanIdx);
        acronym = ibl_id_to_acronym(regionId, regionTable);
        parent = ibl_remap_area(acronym, cfg.area_map);
        if isempty(parent)
            continue;
        end
        if ~ibl_is_allen_consistent_mapping(acronym, parent, regionTable)
            continue;
        end

        st = spikeTimes(spikeClus == cid);
        st = st(st >= spStart & st < spEnd);
        fr = numel(st) / spDur;
        if fr < cfg.min_mean_fr
            continue;
        end

        item.times = st(:)' - spStart;
        item.fr = fr;

        if ~areaBuckets.isKey(parent)
            areaBuckets(parent) = {item};
        else
            tmp = areaBuckets(parent);
            tmp{end+1} = item;
            areaBuckets(parent) = tmp;
        end
    end
end

areas = sort(keys(areaBuckets));
areaStruct = struct('area', {}, 'spike_times', {}, 'firing_rates', {}, 'n_neurons', {});

for k = 1:numel(areas)
    bucket = areaBuckets(areas{k});
    if numel(bucket) < cfg.min_neurons_per_area
        continue;
    end

    fr = cellfun(@(x) x.fr, bucket);
    [~, ord] = sort(fr, 'descend');
    ord = ord(1:min(numel(ord), cfg.max_neurons_per_area));
    bucket = bucket(ord);
    fr = fr(ord);

    if mean(fr) < cfg.min_mean_fr
        continue;
    end

    s.area = string(areas{k});
    s.spike_times = cellfun(@(x) x.times, bucket, 'UniformOutput', false);
    s.firing_rates = fr(:);
    s.n_neurons = numel(bucket);
    areaStruct(end+1) = s; %#ok<AGROW>
end
end


function acronym = ibl_id_to_acronym(regionId, regionTable)
idx = find(regionTable.id == regionId, 1);
if isempty(idx)
    acronym = 'void';
else
    acronym = char(regionTable.acronym(idx));
end
end


function goodMask = ibl_load_good_clusters(metricsFile)
persistent warnedParquet
goodMask = [];

if exist(metricsFile, 'file') ~= 2
    return;
end

if exist('parquetread', 'file') ~= 2
    if isempty(warnedParquet)
        warning(['parquetread is not available; cluster quality filtering is disabled. ', ...
                 'All clusters will be treated as usable.']);
        warnedParquet = true;
    end
    return;
end

T = parquetread(metricsFile);
vars = string(T.Properties.VariableNames);
if any(vars == "label")
    goodMask = T.label == 1;
elseif any(vars == "ks2_label")
    if iscellstr(T.ks2_label) || isstring(T.ks2_label)
        goodMask = string(T.ks2_label) == "good";
    else
        goodMask = T.ks2_label == 1;
    end
end
end
