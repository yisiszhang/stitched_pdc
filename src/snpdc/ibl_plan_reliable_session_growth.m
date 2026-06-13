function plan = ibl_plan_reliable_session_growth(cfg, csdCorrThreshold, varargin)
%IBL_PLAN_RELIABLE_SESSION_GROWTH Greedy reliable session-selection plan.
%
%   plan = ibl_plan_reliable_session_growth(cfg, csdCorrThreshold)
%
% Builds a session inclusion plan before full stitching. A session is
% eligible if it passes basic area/session quality gates. Eligible sessions
% are first restricted to the largest area-overlap component, then a
% CSD-reliable session graph is built. The selected sessions are the
% largest connected component of that reliable graph.
%
% This function only reads cached cross-area CSD summaries. It does not run
% Chronux, matrix completion, Wilson factorization, or PDC.

p = inputParser;
p.addParameter('Band', [1 80], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('Metric', 'coherence_abs', @(x)ischar(x)||isstring(x));
p.addParameter('MinAreas', 2, @(x)isnumeric(x)&&isscalar(x)&&x>=2);
p.addParameter('MinMeanPc1', [], @(x)isnumeric(x)||isempty(x));
p.addParameter('MinOverlapAreas', 2, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('MinOverlapPairs', 1, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('MinNewAreas', 1, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('MaxSessions', inf, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('PreselectLargestOverlapComponent', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('OverlapComponentSelection', 'largest_sessions', @(x)ischar(x)||isstring(x));
p.addParameter('ReliableComponentSelection', 'largest_areas', @(x)ischar(x)||isstring(x));
p.addParameter('RequiredAreas', strings(0,1), @(x)isstring(x)||iscellstr(x)||ischar(x));
p.addParameter('SeedSessionId', "", @(x)isstring(x)||ischar(x));
p.addParameter('SessionFilterFile', "", @(x)isstring(x)||ischar(x));
p.addParameter('SessionIncludeIds', strings(0,1), @(x)isstring(x)||iscellstr(x)||ischar(x));
p.addParameter('SessionExcludeIds', strings(0,1), @(x)isstring(x)||iscellstr(x)||ischar(x));
p.addParameter('OutputFile', "", @(x)isstring(x)||ischar(x));
p.addParameter('MakeFigure', true, @(x)islogical(x)&&isscalar(x));
p.parse(varargin{:});
opt = p.Results;
metric = lower(string(opt.Metric));
if isempty(opt.MinMeanPc1)
    opt.MinMeanPc1 = cfg.min_pc1_explained;
end

sessions = local_load_session_spectra(cfg, opt, metric);
assert(~isempty(sessions), 'No eligible session spectra found in %s.', cfg.cross_spectra_dir);
allEligibleSessions = sessions;
preselectionOverlapGraph = local_session_overlap_graph(allEligibleSessions, opt.MinOverlapAreas);
preselectionComponentTable = local_overlap_component_table(allEligibleSessions, preselectionOverlapGraph);
if opt.PreselectLargestOverlapComponent
    keepIdx = local_choose_overlap_component(allEligibleSessions, preselectionOverlapGraph, opt);
    sessions = sessions(keepIdx);
end
overlapGraph = local_session_overlap_graph(sessions, opt.MinOverlapAreas);
componentTable = local_overlap_component_table(sessions, overlapGraph);

[reliableGraph, pairRows] = local_reliable_session_graph(sessions, csdCorrThreshold, opt);
selectedIdx = local_choose_reliable_component(sessions, reliableGraph, opt);
if isfinite(opt.MaxSessions) && numel(selectedIdx) > opt.MaxSessions
    selectedIdx = selectedIdx(1:opt.MaxSessions);
end
selected = sessions(selectedIdx);
for i = 1:numel(selected)
    selected(i).component_id = 1;
end
traceRows = local_selected_trace_rows(selected);
selectedIds = string({selected.session_id}).';
areaNames = local_union_areas(selected);
[supportCount, pairSupportCount] = local_support_counts(selected, areaNames);
componentIds = [selected.component_id].';
eligibleAreas = local_union_areas(sessions);
excludedQualifiedAreas = setdiff(eligibleAreas, areaNames);

plan.options = opt;
plan.csd_corr_threshold = csdCorrThreshold;
plan.selected_sessions = selected;
plan.selected_session_ids = selectedIds;
plan.selected_table = local_selected_table(selected);
plan.trace_table = local_trace_table(traceRows);
plan.eligible_session_table = local_eligible_table(sessions);
plan.all_eligible_session_table = local_eligible_table(allEligibleSessions);
plan.preselection_overlap_graph = preselectionOverlapGraph;
plan.preselection_component_table = preselectionComponentTable;
plan.overlap_graph = overlapGraph;
plan.overlap_component_table = componentTable;
plan.reliable_graph = reliableGraph;
plan.reliable_pair_table = pairRows;
plan.reliable_component_table = local_overlap_component_table(sessions, reliableGraph);
plan.area_names = areaNames;
plan.eligible_area_names = eligibleAreas;
plan.excluded_qualified_area_names = excludedQualifiedAreas;
plan.component_ids = componentIds;
plan.n_components = numel(unique(componentIds));
plan.area_support_count = supportCount;
plan.pair_support_count = pairSupportCount;
plan.completed_pair_fraction = local_completed_pair_fraction(pairSupportCount);
plan.created_at = string(datetime('now'));

if strlength(string(opt.OutputFile)) > 0
    outFile = string(opt.OutputFile);
    outDir = fileparts(char(outFile));
    if ~isempty(outDir) && exist(outDir, 'dir') ~= 7
        mkdir(outDir);
    end
    save(outFile, 'plan', '-v7.3');
    fprintf('[ibl_plan_reliable_session_growth] saved %s\n', char(outFile));
end

fprintf(['[ibl_plan_reliable_session_growth] selected sessions=%d  areas=%d  ' ...
    'components=%d  completed_pair_fraction=%.3f\n'], ...
    numel(selected), numel(areaNames), plan.n_components, plan.completed_pair_fraction);
if opt.PreselectLargestOverlapComponent && ~isempty(preselectionComponentTable)
    keptHasCA1 = any(plan.eligible_area_names == "CA1");
    allHasCA1 = any(local_union_areas(allEligibleSessions) == "CA1");
    fprintf('[ibl_plan_reliable_session_growth] preselection components=%d  all_has_CA1=%d  kept_has_CA1=%d\n', ...
        height(preselectionComponentTable), allHasCA1, keptHasCA1);
end

if opt.MakeFigure
    local_plot_plan(plan);
end
end


function sessions = local_load_session_spectra(cfg, opt, metric)
files = dir(fullfile(cfg.cross_spectra_dir, '*.mat'));
assert(~isempty(files), 'No cross-spectrum summaries found in %s.', cfg.cross_spectra_dir);
sessionKeep = local_session_keep_set(cfg, opt);
sessions = struct('session_id', {}, 'lab', {}, 'subject', {}, 'date', {}, ...
    'area_names', {}, 'mean_pc1', {}, 'n_areas', {}, 'freqs', {}, ...
    'pair_keys', {}, 'pair_vectors', {}, 'file', {}, 'component_id', {});

fprintf('[ibl_plan_reliable_session_growth] scanning %d CSD files\n', numel(files));
for k = 1:numel(files)
    tmp = load(fullfile(files(k).folder, files(k).name), 'summary');
    summary = tmp.summary;
    sessionId = string(summary.session_id);
    if ~local_keep_session(sessionId, sessionKeep)
        continue;
    end
    if ~isfield(summary, 'cross_spectrum') || ~isfield(summary, 'freqs')
        continue;
    end

    areaScores = local_area_scores(summary, opt.MinMeanPc1);
    keepArea = areaScores >= opt.MinMeanPc1;
    areaNames = string(summary.area_names(keepArea));
    if numel(areaNames) < opt.MinAreas
        continue;
    end
    S = summary.cross_spectrum(keepArea, keepArea, :);
    freqs = summary.freqs(:);
    fMask = freqs >= opt.Band(1) & freqs <= min(opt.Band(2), freqs(end));
    if ~any(fMask)
        continue;
    end

    [areaNames, ord] = sort(areaNames);
    S = S(ord, ord, :);
    [pairKeys, pairVectors] = local_pair_vectors(S, areaNames, fMask, metric);
    if isempty(pairKeys)
        continue;
    end

    parts = split(sessionId, "__");
    rec.session_id = sessionId;
    rec.lab = local_part(parts, 1);
    rec.subject = local_part(parts, 2);
    rec.date = local_part(parts, 3);
    rec.area_names = areaNames(:);
    rec.mean_pc1 = local_nanmean(areaScores(keepArea));
    rec.n_areas = numel(areaNames);
    rec.freqs = freqs(fMask);
    rec.pair_keys = pairKeys;
    rec.pair_vectors = pairVectors;
    rec.file = string(fullfile(files(k).folder, files(k).name));
    rec.component_id = NaN;
    sessions(end+1) = rec; %#ok<AGROW>
end
end


function keepIdx = local_choose_overlap_component(sessions, graph, opt)
mode = lower(string(opt.OverlapComponentSelection));
switch mode
    case {"largest_sessions", "sessions"}
        keepIdx = graph.components{1};
    case {"largest_areas", "areas"}
        nAreas = zeros(numel(graph.components), 1);
        for c = 1:numel(graph.components)
            nAreas(c) = numel(local_union_areas(sessions(graph.components{c})));
        end
        [~, best] = max(nAreas);
        keepIdx = graph.components{best};
    case {"none", "all"}
        keepIdx = 1:numel(sessions);
    otherwise
        error('Unknown OverlapComponentSelection: %s', mode);
end
end


function graph = local_session_overlap_graph(sessions, minOverlapAreas)
n = numel(sessions);
adj = false(n);
overlapCount = zeros(n);
for i = 1:n
    for j = (i+1):n
        nOverlap = numel(intersect(sessions(i).area_names(:), sessions(j).area_names(:)));
        overlapCount(i,j) = nOverlap;
        overlapCount(j,i) = nOverlap;
        if nOverlap >= minOverlapAreas
            adj(i,j) = true;
            adj(j,i) = true;
        end
    end
end
adj(1:n+1:end) = true;

components = local_connected_components(adj);
sessionIds = string({sessions.session_id}).';
componentId = zeros(n,1);
for c = 1:numel(components)
    componentId(components{c}) = c;
end

graph.session_ids = sessionIds;
graph.adjacency = adj;
graph.overlap_count = overlapCount;
graph.components = components;
graph.component_id = componentId;
graph.component_sizes = cellfun(@numel, components(:));
end


function [graph, pairTable] = local_reliable_session_graph(sessions, threshold, opt)
n = numel(sessions);
adj = false(n);
overlapCount = zeros(n);
pairCount = zeros(n);
corrMat = nan(n);
rows = struct('session_i', {}, 'session_j', {}, 'n_overlap_areas', {}, ...
    'n_overlap_pairs', {}, 'median_csd_corr', {}, 'pass', {});

for i = 1:n
    for j = (i+1):n
        [nOverlapAreas, nOverlapPairs, corrVal] = local_pairwise_session_csd_corr(sessions(i), sessions(j));
        overlapCount(i,j) = nOverlapAreas;
        overlapCount(j,i) = nOverlapAreas;
        pairCount(i,j) = nOverlapPairs;
        pairCount(j,i) = nOverlapPairs;
        corrMat(i,j) = corrVal;
        corrMat(j,i) = corrVal;
        pass = nOverlapAreas >= opt.MinOverlapAreas && ...
            nOverlapPairs >= opt.MinOverlapPairs && ...
            isfinite(corrVal) && corrVal >= threshold;
        adj(i,j) = pass;
        adj(j,i) = pass;

        row.session_i = sessions(i).session_id;
        row.session_j = sessions(j).session_id;
        row.n_overlap_areas = nOverlapAreas;
        row.n_overlap_pairs = nOverlapPairs;
        row.median_csd_corr = corrVal;
        row.pass = pass;
        rows(end+1) = row; %#ok<AGROW>
    end
end
adj(1:n+1:end) = true;

components = local_connected_components(adj);
sessionIds = string({sessions.session_id}).';
componentId = zeros(n,1);
for c = 1:numel(components)
    componentId(components{c}) = c;
end

graph.session_ids = sessionIds;
graph.adjacency = adj;
graph.overlap_count = overlapCount;
graph.pair_count = pairCount;
graph.csd_corr = corrMat;
graph.components = components;
graph.component_id = componentId;
graph.component_sizes = cellfun(@numel, components(:));

if isempty(rows)
    pairTable = table();
else
    pairTable = struct2table(rows);
end
end


function selectedIdx = local_choose_reliable_component(sessions, graph, opt)
if isempty(graph.components)
    selectedIdx = [];
    return;
end
mode = lower(string(opt.ReliableComponentSelection));
requiredAreas = string(opt.RequiredAreas(:));
best = 1;
bestScore = [-inf, -inf, -inf];
for c = 1:numel(graph.components)
    idx = graph.components{c};
    areas = local_union_areas(sessions(idx));
    reqHits = numel(intersect(areas(:), requiredAreas));
    switch mode
        case {"largest_sessions", "sessions"}
            score = [reqHits, numel(idx), numel(areas)];
        case {"largest_areas", "areas"}
            score = [reqHits, numel(areas), numel(idx)];
        otherwise
            error('Unknown ReliableComponentSelection: %s', mode);
    end
    if local_lex_gt(score, bestScore)
        bestScore = score;
        best = c;
    end
end
selectedIdx = graph.components{best};
end


function [nOverlapAreas, nOverlapPairs, corrVal] = local_pairwise_session_csd_corr(a, b)
sharedAreas = intersect(a.area_names(:), b.area_names(:));
nOverlapAreas = numel(sharedAreas);
vals = [];
if nOverlapAreas >= 2
    pairKeys = local_area_pair_keys(sharedAreas);
    for p = 1:numel(pairKeys)
        va = local_get_pair_vector(a, pairKeys(p));
        vb = local_get_pair_vector(b, pairKeys(p));
        r = local_vec_corr(va, vb);
        if isfinite(r)
            vals(end+1,1) = r; %#ok<AGROW>
        end
    end
end
nOverlapPairs = numel(vals);
corrVal = local_nanmedian(vals);
end


function rows = local_selected_trace_rows(selected)
rows = struct('step', {}, 'session_id', {}, 'decision', {}, ...
    'component_id', {}, 'n_selected_areas_before', {}, ...
    'n_session_areas', {}, 'n_overlap_areas', {}, 'n_new_areas', {}, ...
    'n_overlap_pairs', {}, 'overlap_csd_corr', {});
areasSoFar = strings(0,1);
for i = 1:numel(selected)
    overlapAreas = intersect(selected(i).area_names(:), areasSoFar(:));
    newAreas = setdiff(selected(i).area_names(:), areasSoFar(:));
    rows(end+1) = local_trace_row(i, selected(i), "reliable_component", 1, ...
        numel(areasSoFar), selected(i).n_areas, numel(overlapAreas), ...
        numel(newAreas), NaN, NaN); %#ok<AGROW>
    areasSoFar = union(areasSoFar, selected(i).area_names(:));
end
end


function T = local_overlap_component_table(sessions, graph)
if isempty(graph.components)
    T = table();
    return;
end
n = numel(graph.components);
component_id = (1:n).';
n_sessions = zeros(n,1);
n_areas = zeros(n,1);
has_CA1 = false(n,1);
session_ids = strings(n,1);
area_names = strings(n,1);
for c = 1:n
    idx = graph.components{c};
    compSessions = sessions(idx);
    areas = local_union_areas(compSessions);
    n_sessions(c) = numel(idx);
    n_areas(c) = numel(areas);
    has_CA1(c) = any(areas == "CA1");
    session_ids(c) = strjoin(string({compSessions.session_id}), ", ");
    area_names(c) = strjoin(areas(:).', ", ");
end
T = table(component_id, n_sessions, n_areas, has_CA1, session_ids, area_names);
end


function components = local_connected_components(adj)
n = size(adj,1);
seen = false(n,1);
components = {};
for i = 1:n
    if seen(i)
        continue;
    end
    queue = i;
    seen(i) = true;
    comp = [];
    while ~isempty(queue)
        cur = queue(1);
        queue(1) = [];
        comp(end+1) = cur; %#ok<AGROW>
        nb = find(adj(cur,:) & ~seen.');
        seen(nb) = true;
        queue = [queue nb]; %#ok<AGROW>
    end
    components{end+1} = comp; %#ok<AGROW>
end
[~, ord] = sort(cellfun(@numel, components), 'descend');
components = components(ord);
end


function row = local_trace_row(step, session, decision, componentId, nBefore, nSessionAreas, nOverlapAreas, nNewAreas, nOverlapPairs, corrVal)
row.step = step;
row.session_id = session.session_id;
row.decision = string(decision);
row.component_id = componentId;
row.n_selected_areas_before = nBefore;
row.n_session_areas = nSessionAreas;
row.n_overlap_areas = nOverlapAreas;
row.n_new_areas = nNewAreas;
row.n_overlap_pairs = nOverlapPairs;
row.overlap_csd_corr = corrVal;
end


function [pairKeys, pairVectors] = local_pair_vectors(S, areaNames, fMask, metric)
pairKeys = strings(0,1);
pairVectors = {};
for i = 1:numel(areaNames)
    for j = (i+1):numel(areaNames)
        v = local_pair_vector(S, i, j, fMask, metric);
        if isempty(v) || all(~isfinite(v))
            continue;
        end
        pairKeys(end+1,1) = string(areaNames(i) + "--" + areaNames(j)); %#ok<AGROW>
        pairVectors{end+1,1} = v(:); %#ok<AGROW>
    end
end
end


function keys = local_area_pair_keys(areaNames)
areaNames = sort(string(areaNames(:)));
keys = strings(0,1);
for i = 1:numel(areaNames)
    for j = (i+1):numel(areaNames)
        keys(end+1,1) = string(areaNames(i) + "--" + areaNames(j)); %#ok<AGROW>
    end
end
end


function v = local_get_pair_vector(session, pairKey)
idx = find(session.pair_keys == pairKey, 1);
if isempty(idx)
    v = [];
else
    v = session.pair_vectors{idx};
end
end


function v = local_pair_vector(S, i, j, fMask, metric)
Sij = squeeze(S(i,j,fMask));
switch metric
    case {"coherence_abs", "coh_abs", "coherence"}
        Sii = squeeze(real(S(i,i,fMask)));
        Sjj = squeeze(real(S(j,j,fMask)));
        denom = sqrt(max(Sii, 0) .* max(Sjj, 0));
        denom(denom <= 0 | ~isfinite(denom)) = NaN;
        v = abs(Sij) ./ denom;
    case {"csd_abs", "abs_csd"}
        v = abs(Sij);
    case {"realimag", "complex"}
        v = [real(Sij(:)); imag(Sij(:))];
    case {"phase"}
        v = angle(Sij);
    otherwise
        error('Unknown metric: %s', metric);
end
v = real(v(:));
end


function tf = local_lex_gt(a, b)
tf = false;
for i = 1:numel(a)
    if a(i) > b(i)
        tf = true;
        return;
    elseif a(i) < b(i)
        return;
    end
end
end


function [supportCount, pairSupportCount] = local_support_counts(selected, areaNames)
n = numel(areaNames);
supportCount = zeros(n,1);
pairSupportCount = zeros(n);
for s = 1:numel(selected)
    [tf, idx] = ismember(selected(s).area_names(:), areaNames);
    idx = idx(tf);
    supportCount(idx) = supportCount(idx) + 1;
    if numel(idx) >= 2
        pairSupportCount(idx, idx) = pairSupportCount(idx, idx) + 1;
    end
end
end


function frac = local_completed_pair_fraction(pairSupportCount)
if isempty(pairSupportCount)
    frac = NaN;
    return;
end
off = triu(~eye(size(pairSupportCount)), 1);
frac = nnz(pairSupportCount(off) == 0) / max(nnz(off), 1);
end


function areas = local_union_areas(selected)
areas = strings(0,1);
for s = 1:numel(selected)
    areas = union(areas, selected(s).area_names(:));
end
areas = sort(areas);
end


function T = local_selected_table(selected)
if isempty(selected)
    T = table();
    return;
end
n = numel(selected);
step = (1:n).';
component_id = zeros(n,1);
session_id = strings(n,1);
lab = strings(n,1);
subject = strings(n,1);
date = strings(n,1);
n_areas = zeros(n,1);
mean_pc1 = nan(n,1);
area_names = strings(n,1);
for i = 1:n
    if isfield(selected, 'component_id')
        component_id(i) = selected(i).component_id;
    else
        component_id(i) = 1;
    end
    session_id(i) = selected(i).session_id;
    lab(i) = selected(i).lab;
    subject(i) = selected(i).subject;
    date(i) = selected(i).date;
    n_areas(i) = selected(i).n_areas;
    mean_pc1(i) = selected(i).mean_pc1;
    area_names(i) = strjoin(selected(i).area_names(:).', ", ");
end
T = table(step, component_id, session_id, lab, subject, date, n_areas, mean_pc1, area_names);
end


function T = local_trace_table(rows)
if isempty(rows)
    T = table();
else
    T = struct2table(rows);
end
end


function T = local_eligible_table(sessions)
if isempty(sessions)
    T = table();
    return;
end
n = numel(sessions);
session_id = strings(n,1);
lab = strings(n,1);
subject = strings(n,1);
n_areas = zeros(n,1);
mean_pc1 = nan(n,1);
area_names = strings(n,1);
for i = 1:n
    session_id(i) = sessions(i).session_id;
    lab(i) = sessions(i).lab;
    subject(i) = sessions(i).subject;
    n_areas(i) = sessions(i).n_areas;
    mean_pc1(i) = sessions(i).mean_pc1;
    area_names(i) = strjoin(sessions(i).area_names(:).', ", ");
end
T = table(session_id, lab, subject, n_areas, mean_pc1, area_names);
T = sortrows(T, {'n_areas', 'mean_pc1'}, {'descend', 'descend'});
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


function scores = local_area_scores(summary, minPc1)
if isfield(summary, 'mean_pc1_explained')
    scores = summary.mean_pc1_explained(:);
else
    scores = nan(numel(summary.area_names), 1);
end
scores(~isfinite(scores)) = -inf;
if isempty(minPc1)
    minPc1 = -inf;
end
end


function r = local_vec_corr(a, b)
a = real(a(:));
b = real(b(:));
good = isfinite(a) & isfinite(b);
if nnz(good) < 3 || std(a(good)) == 0 || std(b(good)) == 0
    r = NaN;
else
    r = corr(a(good), b(good), 'type', 'Pearson');
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


function local_plot_plan(plan)
figure('Color', 'w', 'Position', [100 100 1220 430]);
tiledlayout(1,3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
T = plan.trace_table;
plot(T.step, cumsum(T.n_new_areas), 'o-', 'LineWidth', 1.5);
xlabel('Growth step');
ylabel('Cumulative new areas');
title('Network area growth');
grid on;

nexttile;
bar(plan.area_support_count);
set(gca, 'XTick', 1:numel(plan.area_names), 'XTickLabel', plan.area_names, ...
    'XTickLabelRotation', 90, 'TickLabelInterpreter', 'none');
ylabel('Session support count');
title('Area support');
grid on;

nexttile;
imagesc(plan.pair_support_count);
axis image;
colorbar;
set(gca, 'XTick', 1:numel(plan.area_names), 'XTickLabel', plan.area_names, ...
    'YTick', 1:numel(plan.area_names), 'YTickLabel', plan.area_names, ...
    'XTickLabelRotation', 90, 'TickLabelInterpreter', 'none');
title('Pair support count');
end
