function report = ibl_diagnose_session_neuron_counts(sessionRoot, cfg)
%IBL_DIAGNOSE_SESSION_NEURON_COUNTS Audit why raw clusters differ from PCA counts.
%
%   report = ibl_diagnose_session_neuron_counts(sessionRoot, cfg)
%
% The PCA/CSD pipeline does not use all Kilosort clusters. It keeps only
% good clusters in mapped Allen-consistent areas, active during the passive
% interval, passing the firing-rate threshold, and then caps each area to
% cfg.max_neurons_per_area top firing-rate neurons.

if nargin < 2 || isempty(cfg)
    cfg = ibl_default_config();
end

regionTable = ibl_load_region_table(cfg.cache_root);
interval = ibl_read_spontaneous_interval(fullfile(sessionRoot, 'alf', '_ibl_passivePeriods.intervalsTable.csv'));
assert(~isempty(interval), 'Could not read spontaneous interval for %s', sessionRoot);
spStart = interval(1);
spEnd = interval(2);
spDur = spEnd - spStart;

probeDirs = dir(fullfile(sessionRoot, 'alf', 'probe*', 'pykilosort'));
areaBuckets = containers.Map('KeyType', 'char', 'ValueType', 'any');
probeRows = struct([]);

for p = 1:numel(probeDirs)
    baseDir = fullfile(probeDirs(p).folder, probeDirs(p).name);
    probeName = string(probeDirs(p).folder);
    parts = split(probeName, filesep);
    probeName = parts(end);

    spikesTimesFile = fullfile(baseDir, 'spikes.times.npy');
    spikesClusFile = fullfile(baseDir, 'spikes.clusters.npy');
    clusterChanFile = fullfile(baseDir, 'clusters.channels.npy');
    brainIdsFile = fullfile(baseDir, 'channels.brainLocationIds_ccf_2017.npy');

    row = local_empty_probe_row(probeName);
    if ~all(cellfun(@(f) exist(f, 'file') == 2, ...
            {spikesTimesFile, spikesClusFile, clusterChanFile, brainIdsFile}))
        row.status = "missing_required_files";
        probeRows = local_append_struct(probeRows, row);
        continue;
    end

    spikeTimes = double(ibl_read_npy(spikesTimesFile));
    spikeClus = double(ibl_read_npy(spikesClusFile));
    cluChan = double(ibl_read_npy(clusterChanFile));
    brainIds = double(ibl_read_npy(brainIdsFile));
    nClusters = numel(cluChan);

    goodMask = local_load_good_clusters(fullfile(baseDir, 'clusters.metrics.pqt'), nClusters);
    clusterIds = find(goodMask) - 1;

    inWindow = spikeTimes >= spStart & spikeTimes < spEnd & ...
        spikeClus >= 0 & spikeClus < nClusters;
    if any(inWindow)
        passiveCounts = accumarray(spikeClus(inWindow) + 1, 1, [nClusters 1], @sum, 0);
    else
        passiveCounts = zeros(nClusters, 1);
    end
    firingRates = passiveCounts ./ spDur;

    nMapped = 0;
    nAllenConsistent = 0;
    nFrPass = 0;
    nAreaExcluded = 0;
    nBadChannel = 0;

    for c = 1:numel(clusterIds)
        cid = clusterIds(c);
        chanIdx = cluChan(cid + 1) + 1;
        if chanIdx < 1 || chanIdx > numel(brainIds)
            nBadChannel = nBadChannel + 1;
            continue;
        end

        acronym = local_id_to_acronym(brainIds(chanIdx), regionTable);
        parent = ibl_remap_area(acronym, cfg.area_map);
        if isempty(parent)
            nAreaExcluded = nAreaExcluded + 1;
            continue;
        end
        nMapped = nMapped + 1;

        if ~ibl_is_allen_consistent_mapping(acronym, parent, regionTable)
            continue;
        end
        nAllenConsistent = nAllenConsistent + 1;

        fr = firingRates(cid + 1);
        if fr < cfg.min_mean_fr
            continue;
        end
        nFrPass = nFrPass + 1;

        item.cluster_id = cid;
        item.probe = probeName;
        item.raw_acronym = string(acronym);
        item.fr = fr;
        if ~areaBuckets.isKey(parent)
            areaBuckets(parent) = {item};
        else
            tmp = areaBuckets(parent);
            tmp{end+1} = item;
            areaBuckets(parent) = tmp;
        end
    end

    row.status = "pass";
    row.n_raw_clusters = nClusters;
    row.n_good_clusters = nnz(goodMask);
    row.n_bad_channel = nBadChannel;
    row.n_area_excluded_or_unmapped = nAreaExcluded;
    row.n_mapped = nMapped;
    row.n_allen_consistent = nAllenConsistent;
    row.n_fr_pass = nFrPass;
    probeRows = local_append_struct(probeRows, row);
end

[areaPreCapRows, areaFinalRows] = local_area_rows(areaBuckets, cfg);

report.session_root = string(sessionRoot);
report.spontaneous_interval = interval(:).';
report.spontaneous_duration_s = spDur;
report.config_thresholds = struct( ...
    'min_mean_fr', cfg.min_mean_fr, ...
    'min_neurons_per_area', cfg.min_neurons_per_area, ...
    'max_neurons_per_area', cfg.max_neurons_per_area, ...
    'min_pc1_explained', cfg.min_pc1_explained);
report.probe_table = local_struct_to_table(probeRows);
report.area_pre_cap_table = local_struct_to_table(areaPreCapRows);
report.area_final_table = local_struct_to_table(areaFinalRows);
report.cascade_table = local_cascade_table(report.probe_table, report.area_final_table, cfg);

if ~isempty(report.area_final_table)
    report.n_final_pca_neurons = sum(report.area_final_table.n_neurons_after_cap, 'omitnan');
else
    report.n_final_pca_neurons = 0;
end

fprintf('[ibl_diagnose_session_neuron_counts] raw=%d  good=%d  fr-pass=%d  final-pca=%d\n', ...
    sum(report.probe_table.n_raw_clusters, 'omitnan'), ...
    sum(report.probe_table.n_good_clusters, 'omitnan'), ...
    sum(report.probe_table.n_fr_pass, 'omitnan'), ...
    report.n_final_pca_neurons);
disp(report.cascade_table);
end

function row = local_empty_probe_row(probeName)
row.probe = string(probeName);
row.status = "";
row.n_raw_clusters = 0;
row.n_good_clusters = 0;
row.n_bad_channel = 0;
row.n_area_excluded_or_unmapped = 0;
row.n_mapped = 0;
row.n_allen_consistent = 0;
row.n_fr_pass = 0;
end

function goodMask = local_load_good_clusters(metricsFile, nClusters)
goodMask = true(nClusters, 1);
if exist(metricsFile, 'file') ~= 2 || exist('parquetread', 'file') ~= 2
    return
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
else
    goodMask = true(height(T), 1);
end

goodMask = logical(goodMask(:));
if numel(goodMask) > nClusters
    goodMask = goodMask(1:nClusters);
elseif numel(goodMask) < nClusters
    goodMask(end+1:nClusters) = false;
end
end

function acronym = local_id_to_acronym(regionId, regionTable)
idx = find(regionTable.id == regionId, 1);
if isempty(idx)
    acronym = 'void';
else
    acronym = char(regionTable.acronym(idx));
end
end

function [preRows, finalRows] = local_area_rows(areaBuckets, cfg)
areas = sort(keys(areaBuckets));
preRows = struct([]);
finalRows = struct([]);

for k = 1:numel(areas)
    area = string(areas{k});
    bucket = areaBuckets(areas{k});
    fr = cellfun(@(x) x.fr, bucket);
    probes = strings(numel(bucket), 1);
    rawAcronyms = strings(numel(bucket), 1);
    for i = 1:numel(bucket)
        probes(i) = bucket{i}.probe;
        rawAcronyms(i) = bucket{i}.raw_acronym;
    end

    pre.area = area;
    pre.n_neurons_before_area_min = numel(bucket);
    pre.mean_fr_before_cap = mean(fr, 'omitnan');
    pre.max_fr = max(fr, [], 'omitnan');
    pre.probes = strjoin(unique(probes, 'stable').', "; ");
    pre.raw_acronyms = strjoin(unique(rawAcronyms, 'stable').', "; ");
    preRows = local_append_struct(preRows, pre);

    if numel(bucket) < cfg.min_neurons_per_area
        final = local_final_area_row(area, numel(bucket), 0, mean(fr, 'omitnan'), ...
            "too_few_neurons_per_area");
        finalRows = local_append_struct(finalRows, final);
        continue;
    end

    [~, ord] = sort(fr, 'descend');
    ord = ord(1:min(numel(ord), cfg.max_neurons_per_area));
    frKept = fr(ord);

    if mean(frKept, 'omitnan') < cfg.min_mean_fr
        final = local_final_area_row(area, numel(bucket), 0, mean(frKept, 'omitnan'), ...
            "low_mean_area_fr_after_cap");
    else
        final = local_final_area_row(area, numel(bucket), numel(frKept), ...
            mean(frKept, 'omitnan'), "used_for_pca");
    end
    finalRows = local_append_struct(finalRows, final);
end
end

function row = local_final_area_row(area, nBefore, nAfter, meanFr, status)
row.area = area;
row.n_neurons_before_cap = nBefore;
row.n_neurons_after_cap = nAfter;
row.mean_fr_after_cap = meanFr;
row.status = string(status);
end

function rows = local_append_struct(rows, row)
if isempty(rows)
    rows = row;
else
    rows(end+1) = row; %#ok<AGROW>
end
end

function T = local_struct_to_table(rows)
if isempty(rows)
    T = table();
else
    T = struct2table(rows);
end
end

function T = local_cascade_table(probeTable, areaFinalTable, cfg)
stage = strings(0,1);
n_neurons = zeros(0,1);
removed_from_previous = zeros(0,1);
note = strings(0,1);

raw = local_sum_table_var(probeTable, 'n_raw_clusters');
good = local_sum_table_var(probeTable, 'n_good_clusters');
validChannelGood = good - local_sum_table_var(probeTable, 'n_bad_channel');
mapped = local_sum_table_var(probeTable, 'n_mapped');
allenConsistent = local_sum_table_var(probeTable, 'n_allen_consistent');
frPass = local_sum_table_var(probeTable, 'n_fr_pass');

if isempty(areaFinalTable)
    areaMinPass = 0;
    finalCount = 0;
else
    usable = areaFinalTable.status == "used_for_pca";
    areaMinPass = sum(areaFinalTable.n_neurons_before_cap(usable), 'omitnan');
    finalCount = sum(areaFinalTable.n_neurons_after_cap(usable), 'omitnan');
end

[stage, n_neurons, removed_from_previous, note] = local_add_cascade_row( ...
    stage, n_neurons, removed_from_previous, note, ...
    "raw_clusters", raw, nan, ...
    "all clusters in clusters.channels.npy across available probes");
[stage, n_neurons, removed_from_previous, note] = local_add_cascade_row( ...
    stage, n_neurons, removed_from_previous, note, ...
    "good_clusters", good, raw, ...
    "keeps clusters.metrics label==1 or ks2_label==good");
[stage, n_neurons, removed_from_previous, note] = local_add_cascade_row( ...
    stage, n_neurons, removed_from_previous, note, ...
    "good_clusters_with_valid_channel", validChannelGood, good, ...
    "drops clusters whose assigned channel is outside channels.brainLocationIds");
[stage, n_neurons, removed_from_previous, note] = local_add_cascade_row( ...
    stage, n_neurons, removed_from_previous, note, ...
    "mapped_standardized_area", mapped, validChannelGood, ...
    "drops excluded/unmapped Allen acronyms such as fiber tracts or void");
[stage, n_neurons, removed_from_previous, note] = local_add_cascade_row( ...
    stage, n_neurons, removed_from_previous, note, ...
    "allen_consistent_mapping", allenConsistent, mapped, ...
    "keeps mappings consistent with the Allen hierarchy");
[stage, n_neurons, removed_from_previous, note] = local_add_cascade_row( ...
    stage, n_neurons, removed_from_previous, note, ...
    "spontaneous_fr_pass", frPass, allenConsistent, ...
    sprintf("keeps neurons with passive-period firing rate >= %.3g Hz", cfg.min_mean_fr));
[stage, n_neurons, removed_from_previous, note] = local_add_cascade_row( ...
    stage, n_neurons, removed_from_previous, note, ...
    "area_min_pass_before_cap", areaMinPass, frPass, ...
    sprintf("keeps areas with at least %d retained neurons", cfg.min_neurons_per_area));
[stage, n_neurons, removed_from_previous, note] = local_add_cascade_row( ...
    stage, n_neurons, removed_from_previous, note, ...
    "final_pca_neurons_after_area_cap", finalCount, areaMinPass, ...
    sprintf("caps each area to top %d firing-rate neurons", cfg.max_neurons_per_area));

T = table(stage, n_neurons, removed_from_previous, note);
end

function total = local_sum_table_var(T, varName)
if isempty(T) || ~ismember(varName, T.Properties.VariableNames)
    total = 0;
else
    total = sum(T.(varName), 'omitnan');
end
end

function [stage, n_neurons, removed_from_previous, note] = local_add_cascade_row( ...
    stage, n_neurons, removed_from_previous, note, stageName, nValue, previousValue, rowNote)
stage(end+1,1) = string(stageName);
n_neurons(end+1,1) = nValue;
if isnan(previousValue)
    removed_from_previous(end+1,1) = nan;
else
    removed_from_previous(end+1,1) = previousValue - nValue;
end
note(end+1,1) = string(rowNote);
end
