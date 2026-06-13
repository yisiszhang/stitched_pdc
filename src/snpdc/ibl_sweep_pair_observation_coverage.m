function sweep = ibl_sweep_pair_observation_coverage(plan, thresholds, varargin)
%IBL_SWEEP_PAIR_OBSERVATION_COVERAGE Coverage under observed-pair constraints.
%
%   sweep = ibl_sweep_pair_observation_coverage(plan, [0.3 0.4 0.5 0.6])
%
% Uses the reliable component already selected in PLAN. For each target
% observed-pair fraction, it greedily prunes sessions while preserving
% reliable-graph connectivity, aiming to retain the most brain areas that
% satisfy the requested pair-observation density.

p = inputParser;
p.addParameter('OutputFile', "", @(x)isstring(x)||ischar(x));
p.addParameter('Verbose', true, @(x)islogical(x)&&isscalar(x));
p.parse(varargin{:});
opt = p.Results;

if ischar(plan) || isstring(plan)
    tmp = load(plan, 'plan');
    plan = tmp.plan;
end
if nargin < 2 || isempty(thresholds)
    thresholds = [0.3 0.4 0.5 0.6];
end
thresholds = thresholds(:);

sessions = plan.selected_sessions(:);
assert(~isempty(sessions), 'plan.selected_sessions is empty.');

adj = local_selected_adjacency(plan, sessions);
rows = struct([]);
solutions = struct('target_observed_pair_fraction', {}, 'session_ids', {}, 'area_names', {});

for t = 1:numel(thresholds)
    tau = thresholds(t);
    selectedIdx = local_prune_to_threshold(sessions, adj, tau);
    stats = local_subset_stats(sessions(selectedIdx));

    row.target_observed_pair_fraction = tau;
    row.achieved_observed_pair_fraction = stats.observed_pair_fraction;
    row.missing_pair_fraction = 1 - stats.observed_pair_fraction;
    row.n_sessions = numel(selectedIdx);
    row.n_areas = numel(stats.area_names);
    row.n_total_pairs = stats.n_total_pairs;
    row.n_observed_pairs = stats.n_observed_pairs;
    row.n_missing_pairs = stats.n_total_pairs - stats.n_observed_pairs;
    row.mean_pair_support_observed = stats.mean_pair_support_observed;
    row.median_pair_support_observed = stats.median_pair_support_observed;
    row.max_pair_support_observed = stats.max_pair_support_observed;
    row.connected = local_is_connected(adj(selectedIdx, selectedIdx));
    row.session_ids = strjoin(string({sessions(selectedIdx).session_id}).', "; ");
    row.area_names = strjoin(stats.area_names(:).', "; ");
    rows = local_append_struct(rows, row);

    sol.target_observed_pair_fraction = tau;
    sol.session_ids = string({sessions(selectedIdx).session_id}).';
    sol.area_names = stats.area_names(:);
    solutions(end+1) = sol; %#ok<AGROW>
end

sweep.table = struct2table(rows);
sweep.solutions = solutions;
sweep.thresholds = thresholds;
sweep.source = "plan.selected_sessions";
sweep.created_at = string(datetime('now'));

if strlength(string(opt.OutputFile)) > 0
    outFile = string(opt.OutputFile);
    outDir = fileparts(char(outFile));
    if ~isempty(outDir) && exist(outDir, 'dir') ~= 7
        mkdir(outDir);
    end
    save(outFile, 'sweep', '-v7.3');
    writetable(sweep.table, replace(outFile, ".mat", ".csv"));
end

if opt.Verbose
    disp(sweep.table(:, {'target_observed_pair_fraction', ...
        'achieved_observed_pair_fraction', 'missing_pair_fraction', ...
        'n_sessions', 'n_areas', 'n_observed_pairs', 'n_total_pairs', ...
        'median_pair_support_observed'}));
end
end

function selectedIdx = local_prune_to_threshold(sessions, adj, tau)
selectedIdx = (1:numel(sessions)).';
stats = local_subset_stats(sessions(selectedIdx));
if stats.observed_pair_fraction >= tau
    return;
end

while numel(selectedIdx) > 1
    bestIdx = NaN;
    bestScore = [-inf, -inf, -inf, -inf];
    for k = 1:numel(selectedIdx)
        candIdx = selectedIdx;
        candIdx(k) = [];
        if isempty(candIdx) || ~local_is_connected(adj(candIdx, candIdx))
            continue;
        end
        candStats = local_subset_stats(sessions(candIdx));
        score = [ ...
            candStats.observed_pair_fraction, ...
            numel(candStats.area_names), ...
            candStats.n_observed_pairs, ...
            candStats.median_pair_support_observed];
        if local_lex_gt(score, bestScore)
            bestScore = score;
            bestIdx = k;
        end
    end
    if isnan(bestIdx)
        break;
    end
    selectedIdx(bestIdx) = [];
    stats = local_subset_stats(sessions(selectedIdx));
    if stats.observed_pair_fraction >= tau
        break;
    end
end
end

function stats = local_subset_stats(sessions)
areaNames = local_union_areas(sessions);
n = numel(areaNames);
pairSupport = zeros(n);
for s = 1:numel(sessions)
    [tf, idx] = ismember(string(sessions(s).area_names(:)), areaNames);
    idx = idx(tf);
    if numel(idx) >= 2
        pairSupport(idx, idx) = pairSupport(idx, idx) + 1;
    end
end
off = triu(~eye(n), 1);
supportVals = pairSupport(off);
observedVals = supportVals(supportVals > 0);
stats.area_names = areaNames(:);
stats.pair_support = pairSupport;
stats.n_total_pairs = nnz(off);
stats.n_observed_pairs = nnz(supportVals > 0);
stats.observed_pair_fraction = stats.n_observed_pairs / max(stats.n_total_pairs, 1);
if isempty(observedVals)
    stats.mean_pair_support_observed = NaN;
    stats.median_pair_support_observed = NaN;
    stats.max_pair_support_observed = NaN;
else
    stats.mean_pair_support_observed = mean(observedVals, 'omitnan');
    stats.median_pair_support_observed = median(observedVals, 'omitnan');
    stats.max_pair_support_observed = max(observedVals);
end
end

function adj = local_selected_adjacency(plan, sessions)
n = numel(sessions);
adj = true(n);
if ~isfield(plan, 'reliable_graph') || ~isfield(plan.reliable_graph, 'adjacency')
    return;
end
allIds = string(plan.reliable_graph.session_ids(:));
selIds = string({sessions.session_id}).';
[tf, idx] = ismember(selIds, allIds);
if all(tf)
    adj = logical(plan.reliable_graph.adjacency(idx, idx));
else
    adj = true(n);
end
adj(1:n+1:end) = true;
end

function tf = local_is_connected(adj)
n = size(adj, 1);
if n <= 1
    tf = true;
    return;
end
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

function areas = local_union_areas(sessions)
areas = strings(0,1);
for s = 1:numel(sessions)
    areas = union(areas, string(sessions(s).area_names(:)));
end
areas = sort(areas);
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

function rows = local_append_struct(rows, row)
if isempty(rows)
    rows = row;
else
    rows(end+1) = row; %#ok<AGROW>
end
end
