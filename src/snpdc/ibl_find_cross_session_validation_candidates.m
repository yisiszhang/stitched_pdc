function out = ibl_find_cross_session_validation_candidates(qc, cfg, varargin)
%IBL_FIND_CROSS_SESSION_VALIDATION_CANDIDATES Find source/target validation sets.
%
%   out = ibl_find_cross_session_validation_candidates(qc, cfg)
%
% A target session is held out as the full simultaneous control. Source
% sessions are selected from other animals and stitched to cover as much of
% the target area set as possible. This function only reads cached area-PCA
% summaries and does not compute spectra.

p = inputParser;
p.addParameter('MinTargetAreas', 8, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('MaxTargetAreas', 14, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('MinSourceAreas', 3, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('MinDuration', 600, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('MinMeanPc1', [], @(x)isnumeric(x)||isempty(x));
p.addParameter('MaxSources', 4, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('MinCoverage', 0.7, @(x)isnumeric(x)&&isscalar(x)&&x>0&&x<=1);
p.addParameter('MinAnchorAreas', 2, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('MinAnchorFraction', 0, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('MinObservedPairFraction', 0, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('MaxCompletedPairFraction', 1, @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('RequireConnectedSources', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('MinSourceOverlapAreas', 1, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('RequireDifferentSubject', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('RequireDifferentLab', false, @(x)islogical(x)&&isscalar(x));
p.addParameter('SessionFilterFile', "", @(x)isstring(x)||ischar(x));
p.addParameter('SessionIncludeIds', strings(0,1), @(x)isstring(x)||iscellstr(x)||ischar(x));
p.addParameter('SessionExcludeIds', strings(0,1), @(x)isstring(x)||iscellstr(x)||ischar(x));
p.addParameter('TopN', 30, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('Verbose', true, @(x)islogical(x)&&isscalar(x));
p.parse(varargin{:});
opt = p.Results;

if nargin < 1 || isempty(qc)
    tmp = load(cfg.pca_qc_file, 'qc');
    qc = tmp.qc;
end
if isempty(opt.MinMeanPc1)
    opt.MinMeanPc1 = cfg.min_pc1_explained;
end

sessionInfo = local_load_session_info(qc, cfg, opt);
if isempty(sessionInfo)
    out = struct('candidates', table(), 'sets', struct([]), 'session_info', table(), 'options', opt);
    return;
end

sets = struct([]);
rows = struct('target_session_id', {}, 'target_lab', {}, 'target_subject', {}, ...
    'target_duration_s', {}, 'n_target_areas', {}, 'n_sources', {}, ...
    'source_session_ids', {}, 'source_labs', {}, 'source_subjects', {}, ...
    'n_covered_target_areas', {}, 'coverage', {}, 'n_anchor_areas', {}, ...
    'anchor_fraction', {}, 'n_total_pairs', {}, 'n_observed_pairs', {}, ...
    'n_completed_pairs', {}, 'observed_pair_fraction', {}, ...
    'completed_pair_fraction', {}, 'source_graph_connected', {}, ...
    'min_source_overlap', {}, 'median_source_overlap', {}, ...
    'mean_source_jaccard', {}, 'source_diversity_subjects', {}, ...
    'source_diversity_labs', {}, 'score', {}, 'target_areas', {}, ...
    'covered_areas', {}, 'anchor_areas', {}, 'set_index', {});

for t = 1:numel(sessionInfo)
    target = sessionInfo(t);
    if target.duration_s < opt.MinDuration || numel(target.area_names) < opt.MinTargetAreas
        continue;
    end

    targetAreas = target.area_names(:);
    targetScores = target.area_scores(:);
    [targetScores, ord] = sort(targetScores, 'descend'); %#ok<ASGLU>
    targetAreas = targetAreas(ord);
    if isfinite(opt.MaxTargetAreas)
        keep = 1:min(numel(targetAreas), opt.MaxTargetAreas);
        targetAreas = sort(targetAreas(keep));
    else
        targetAreas = sort(targetAreas);
    end

    sources = local_candidate_sources(sessionInfo, target, targetAreas, opt);
    if isempty(sources)
        continue;
    end

    selected = local_greedy_sources(sources, targetAreas, opt);
    if isempty(selected)
        continue;
    end

    sourceAreas = strings(0,1);
    sourceIds = strings(numel(selected),1);
    sourceLabs = strings(numel(selected),1);
    sourceSubjects = strings(numel(selected),1);
    recset = cell(1, numel(selected));
    for s = 1:numel(selected)
        src = selected(s);
        sourceIds(s) = src.session_id;
        sourceLabs(s) = src.lab;
        sourceSubjects(s) = src.subject;
        srcAreas = intersect(src.area_names(:), targetAreas, 'stable');
        sourceAreas = union(sourceAreas, srcAreas);
        [~, idx] = ismember(srcAreas, targetAreas);
        recset{s} = idx(:).';
    end
    coveredAreas = intersect(targetAreas, sourceAreas);
    coverage = numel(coveredAreas) / numel(targetAreas);
    countMat = meacount_mat(recset, numel(targetAreas));
    observedPair = countMat > 0;
    off = ~eye(numel(targetAreas));
    totalPairs = nnz(triu(off, 1));
    nObservedPairs = nnz(triu(observedPair & off, 1));
    nCompletedPairs = nnz(triu(~observedPair & off, 1));
    anchorAreas = targetAreas(sum(countMat > 0, 2) >= 2);
    nAnchor = numel(anchorAreas);
    anchorFraction = nAnchor / max(numel(coveredAreas), 1);
    observedPairFraction = nObservedPairs / max(totalPairs, 1);
    completedPairFraction = nCompletedPairs / max(totalPairs, 1);
    overlapStats = local_source_overlap_stats(recset, opt.MinSourceOverlapAreas);

    if coverage < opt.MinCoverage || ...
            nAnchor < opt.MinAnchorAreas || ...
            anchorFraction < opt.MinAnchorFraction || ...
            observedPairFraction < opt.MinObservedPairFraction || ...
            completedPairFraction > opt.MaxCompletedPairFraction || ...
            nCompletedPairs < 1 || ...
            (opt.RequireConnectedSources && ~overlapStats.connected)
        continue;
    end

    score = 100 * coverage + 10 * observedPairFraction + 5 * nAnchor + nCompletedPairs + ...
        2 * numel(unique(sourceSubjects)) + numel(unique(sourceLabs));

    setRec.target = target;
    setRec.target_area_names = targetAreas;
    setRec.sources = selected;
    setRec.source_session_ids = sourceIds;
    setRec.source_area_names = arrayfun(@(x) intersect(x.area_names(:), targetAreas, 'stable'), ...
        selected, 'UniformOutput', false);
    setRec.covered_area_names = coveredAreas;
    setRec.anchor_area_names = anchorAreas;
    setRec.count_mat = countMat;
    setRec.observed_pair = observedPair;
    setRec.coverage = coverage;
    setRec.anchor_fraction = anchorFraction;
    setRec.observed_pair_fraction = observedPairFraction;
    setRec.completed_pair_fraction = completedPairFraction;
    setRec.source_overlap = overlapStats;
    setRec.n_completed_pairs = nCompletedPairs;
    sets = [sets; setRec]; %#ok<AGROW>
    setIndex = numel(sets);

    row.target_session_id = target.session_id;
    row.target_lab = target.lab;
    row.target_subject = target.subject;
    row.target_duration_s = target.duration_s;
    row.n_target_areas = numel(targetAreas);
    row.n_sources = numel(selected);
    row.source_session_ids = strjoin(sourceIds.', ", ");
    row.source_labs = strjoin(unique(sourceLabs).', ", ");
    row.source_subjects = strjoin(unique(sourceSubjects).', ", ");
    row.n_covered_target_areas = numel(coveredAreas);
    row.coverage = coverage;
    row.n_anchor_areas = nAnchor;
    row.anchor_fraction = anchorFraction;
    row.n_total_pairs = totalPairs;
    row.n_observed_pairs = nObservedPairs;
    row.n_completed_pairs = nCompletedPairs;
    row.observed_pair_fraction = observedPairFraction;
    row.completed_pair_fraction = completedPairFraction;
    row.source_graph_connected = overlapStats.connected;
    row.min_source_overlap = overlapStats.min_overlap;
    row.median_source_overlap = overlapStats.median_overlap;
    row.mean_source_jaccard = overlapStats.mean_jaccard;
    row.source_diversity_subjects = numel(unique(sourceSubjects));
    row.source_diversity_labs = numel(unique(sourceLabs));
    row.score = score;
    row.target_areas = strjoin(targetAreas.', ", ");
    row.covered_areas = strjoin(coveredAreas.', ", ");
    row.anchor_areas = strjoin(anchorAreas.', ", ");
    row.set_index = setIndex;
    rows(end+1) = row; %#ok<AGROW>
end

if isempty(rows)
    candidateTable = table();
else
    candidateTable = struct2table(rows);
    candidateTable = sortrows(candidateTable, ...
        {'score', 'coverage', 'observed_pair_fraction', 'n_anchor_areas', 'n_completed_pairs'}, ...
        {'descend', 'descend', 'descend', 'descend', 'descend'});
    sets = sets(candidateTable.set_index);
    if isfinite(opt.TopN) && height(candidateTable) > opt.TopN
        candidateTable = candidateTable(1:opt.TopN, :);
        sets = sets(1:height(candidateTable));
    end
    candidateTable.set_index = (1:height(candidateTable)).';
end

out.candidates = candidateTable;
out.sets = sets;
out.session_info = local_session_info_table(sessionInfo);
out.options = opt;

if opt.Verbose
    fprintf('[ibl_find_cross_session_validation_candidates] found %d candidate source/target sets\n', ...
        height(candidateTable));
    if ~isempty(candidateTable)
        disp(candidateTable(:, {'target_session_id', 'n_target_areas', 'n_sources', ...
            'coverage', 'n_anchor_areas', 'observed_pair_fraction', ...
            'completed_pair_fraction', 'source_graph_connected', 'score'}));
    end
end
end


function sessionInfo = local_load_session_info(qc, cfg, opt)
sessions = qc.qualifying_sessions;
sessionKeep = local_session_keep_set(cfg, opt);
sessionInfo = struct('session_id', {}, 'lab', {}, 'subject', {}, 'date', {}, ...
    'number', {}, 'duration_s', {}, 'area_names', {}, 'area_scores', {}, ...
    'mean_pc1', {}, 'pca_file', {});

for k = 1:numel(sessions)
    sessionId = string(sessions(k).session_id);
    if ~local_keep_session(sessionId, sessionKeep)
        continue;
    end
    pcaFile = fullfile(cfg.area_pca_dir, char(sessionId + ".mat"));
    if exist(pcaFile, 'file') ~= 2
        continue;
    end
    tmp = load(pcaFile, 'summary');
    summary = tmp.summary;
    scoresAll = summary.mean_pc1_explained(:);
    keep = summary.pass_qc_area(:) & ...
        ismember(string(summary.area_names(:)), qc.largest_component(:)) & ...
        isfinite(scoresAll) & scoresAll >= opt.MinMeanPc1;
    areaNames = string(summary.area_names(keep));
    areaScores = scoresAll(keep);
    if summary.sp_dur < opt.MinDuration || numel(areaNames) < min(opt.MinSourceAreas, opt.MinTargetAreas)
        continue;
    end
    [areaNames, ord] = sort(areaNames);
    areaScores = areaScores(ord);

    parts = split(sessionId, "__");
    rec.session_id = sessionId;
    rec.lab = local_part(parts, 1);
    rec.subject = local_part(parts, 2);
    rec.date = local_part(parts, 3);
    rec.number = local_part(parts, 4);
    rec.duration_s = summary.sp_dur;
    rec.area_names = areaNames(:);
    rec.area_scores = areaScores(:);
    rec.mean_pc1 = mean(areaScores, 'omitnan');
    rec.pca_file = string(pcaFile);
    sessionInfo(end+1) = rec; %#ok<AGROW>
end
end


function sessionKeep = local_session_keep_set(cfg, opt)
sessionKeep.include = strings(0,1);
sessionKeep.exclude = strings(0,1);

if isfield(cfg, 'session_include_ids') && ~isempty(cfg.session_include_ids)
    sessionKeep.include = string(cfg.session_include_ids(:));
end
if ~isempty(opt.SessionIncludeIds)
    optInclude = string(opt.SessionIncludeIds(:));
    if isempty(sessionKeep.include)
        sessionKeep.include = optInclude;
    else
        sessionKeep.include = intersect(sessionKeep.include, optInclude, 'stable');
    end
end
if isfield(cfg, 'session_exclude_ids') && ~isempty(cfg.session_exclude_ids)
    sessionKeep.exclude = string(cfg.session_exclude_ids(:));
end
if ~isempty(opt.SessionExcludeIds)
    sessionKeep.exclude = unique([sessionKeep.exclude; string(opt.SessionExcludeIds(:))], 'stable');
end

filterFile = string(opt.SessionFilterFile);
if strlength(filterFile) == 0 && isfield(cfg, 'session_filter_file')
    filterFile = string(cfg.session_filter_file);
end
if strlength(filterFile) > 0
    assert(exist(filterFile, 'file') == 2, 'SessionFilterFile not found: %s', char(filterFile));
    tmp = load(filterFile, 'filt');
    assert(isfield(tmp, 'filt') && isfield(tmp.filt, 'kept_session_ids'), ...
        'SessionFilterFile must contain filt.kept_session_ids.');
    filterInclude = string(tmp.filt.kept_session_ids(:));
    if isempty(sessionKeep.include)
        sessionKeep.include = filterInclude;
    else
        sessionKeep.include = intersect(sessionKeep.include, filterInclude, 'stable');
    end
end
end


function tf = local_keep_session(sessionId, sessionKeep)
sessionId = string(sessionId);
if ~isempty(sessionKeep.include) && ~any(sessionKeep.include == sessionId)
    tf = false;
    return;
end
if ~isempty(sessionKeep.exclude) && any(sessionKeep.exclude == sessionId)
    tf = false;
    return;
end
tf = true;
end


function sources = local_candidate_sources(sessionInfo, target, targetAreas, opt)
sources = struct([]);
for i = 1:numel(sessionInfo)
    src = sessionInfo(i);
    if src.session_id == target.session_id
        continue;
    end
    if opt.RequireDifferentSubject && src.subject == target.subject
        continue;
    end
    if opt.RequireDifferentLab && src.lab == target.lab
        continue;
    end
    shared = intersect(src.area_names(:), targetAreas);
    if numel(shared) < opt.MinSourceAreas
        continue;
    end
    src.shared_target_areas = shared(:);
    src.n_shared_target_areas = numel(shared);
    sources = [sources; src]; %#ok<AGROW>
end
end


function T = local_session_info_table(sessionInfo)
if isempty(sessionInfo)
    T = table();
    return;
end
n = numel(sessionInfo);
session_id = strings(n,1);
lab = strings(n,1);
subject = strings(n,1);
date = strings(n,1);
number = strings(n,1);
duration_s = zeros(n,1);
n_areas = zeros(n,1);
mean_pc1 = nan(n,1);
area_names = strings(n,1);
pca_file = strings(n,1);
for i = 1:n
    session_id(i) = sessionInfo(i).session_id;
    lab(i) = sessionInfo(i).lab;
    subject(i) = sessionInfo(i).subject;
    date(i) = sessionInfo(i).date;
    number(i) = sessionInfo(i).number;
    duration_s(i) = sessionInfo(i).duration_s;
    n_areas(i) = numel(sessionInfo(i).area_names);
    mean_pc1(i) = sessionInfo(i).mean_pc1;
    area_names(i) = strjoin(sessionInfo(i).area_names(:).', ", ");
    pca_file(i) = sessionInfo(i).pca_file;
end
T = table(session_id, lab, subject, date, number, duration_s, n_areas, ...
    mean_pc1, area_names, pca_file);
end


function selected = local_greedy_sources(sources, targetAreas, opt)
selected = struct([]);
covered = strings(0,1);
for step = 1:opt.MaxSources
    bestIdx = 0;
    bestScore = -inf;
    for i = 1:numel(sources)
        if any(arrayfun(@(x) x.session_id == sources(i).session_id, selected))
            continue;
        end
        srcAreas = intersect(sources(i).area_names(:), targetAreas);
        newAreas = setdiff(srcAreas, covered);
        anchorBonus = numel(intersect(srcAreas, covered));
        diversityBonus = 0;
        if isempty(selected) || ~any(arrayfun(@(x) x.subject == sources(i).subject, selected))
            diversityBonus = diversityBonus + 1;
        end
        if isempty(selected) || ~any(arrayfun(@(x) x.lab == sources(i).lab, selected))
            diversityBonus = diversityBonus + 0.5;
        end
        score = 10 * numel(newAreas) + 2 * anchorBonus + diversityBonus + ...
            0.01 * sources(i).duration_s + sources(i).mean_pc1;
        if score > bestScore
            bestScore = score;
            bestIdx = i;
        end
    end
    if bestIdx == 0
        break;
    end
    selected = [selected; sources(bestIdx)]; %#ok<AGROW>
    covered = union(covered, intersect(sources(bestIdx).area_names(:), targetAreas));
    if numel(covered) / numel(targetAreas) >= opt.MinCoverage && numel(selected) >= 2
        break;
    end
end
end


function stats = local_source_overlap_stats(recset, minOverlapAreas)
n = numel(recset);
overlap = zeros(n);
jaccard = nan(n);
for i = 1:n
    ai = unique(recset{i});
    for j = (i+1):n
        aj = unique(recset{j});
        inter = numel(intersect(ai, aj));
        uni = numel(union(ai, aj));
        overlap(i,j) = inter;
        overlap(j,i) = inter;
        if uni > 0
            jaccard(i,j) = inter / uni;
            jaccard(j,i) = jaccard(i,j);
        end
    end
end

if n <= 1
    connected = true;
else
    adj = overlap >= minOverlapAreas;
    adj(1:n+1:end) = true;
    connected = local_is_connected(adj);
end

vals = overlap(triu(true(n), 1));
jvals = jaccard(triu(true(n), 1));
stats.connected = connected;
stats.overlap = overlap;
stats.jaccard = jaccard;
stats.min_overlap = local_nanmin(vals);
stats.median_overlap = local_nanmedian(vals);
stats.mean_jaccard = local_nanmean(jvals);
end


function tf = local_is_connected(adj)
n = size(adj, 1);
seen = false(n,1);
queue = 1;
seen(1) = true;
while ~isempty(queue)
    cur = queue(1);
    queue(1) = [];
    nb = find(adj(cur,:) & ~seen.');
    seen(nb) = true;
    queue = [queue nb]; %#ok<AGROW>
end
tf = all(seen);
end


function y = local_nanmin(x)
x = x(isfinite(x));
if isempty(x)
    y = NaN;
else
    y = min(x);
end
end


function y = local_nanmedian(x)
x = x(isfinite(x));
if isempty(x)
    y = NaN;
else
    y = median(x);
end
end


function y = local_nanmean(x)
x = x(isfinite(x));
if isempty(x)
    y = NaN;
else
    y = mean(x);
end
end


function out = local_part(parts, idx)
if numel(parts) >= idx
    out = string(parts(idx));
else
    out = "";
end
end
