%SUMMARIZE_IBL_METHOD_SESSIONS Summarize sessions/areas used in method checks.
%
% Run from the repository root:
%   cd('/Users/metis/Documents/MATLAB/stitch_causality/icml/snpdc_package/snpdc_icml')
%   run('scripts/summarize_ibl_method_sessions.m')
%
% Outputs:
%   ibl_output/method_verification/session_summary/session_summary.csv
%   ibl_output/method_verification/session_summary/session_area_summary.csv
%   ibl_output/method_verification/session_summary/session_method_roles.csv
%   ibl_output/method_verification/session_summary/session_method_summary.mat

clearvars;
repoRoot = pwd;
setup;
cfg = ibl_default_config();

outDir = fullfile(cfg.output_root, 'method_verification', 'session_summary');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

withinFile = fullfile(cfg.output_root, 'validation', ...
    'within_session_subset_stitching_coherence_pdc_20repeats.mat');
crossFile = fullfile(cfg.output_root, 'validation', ...
    'cross_session_filtered_overlap_pdc.mat');
planFile = fullfile(cfg.output_root, 'reliable_session_component_plan.mat');
if ~exist(planFile, 'file')
    planFile = fullfile(cfg.output_root, 'reliable_session_growth_plan.mat');
end
finalFile = fullfile(cfg.output_root, 'reliable_component_result', ...
    'stitched_pdc_reliable_component.mat');

roleRows = local_collect_method_roles(withinFile, crossFile, planFile, finalFile);
roleTable = local_roles_to_table(roleRows);
roleTable = unique(roleTable, 'rows', 'stable');

finalSessionIds = local_final_session_ids(planFile, finalFile);
finalAreaNames = local_final_area_names(finalFile);

sessionIds = unique(roleTable.session_id, 'stable');
[sessionRows, areaRows] = local_build_summary_rows(sessionIds, roleTable, ...
    finalSessionIds, finalAreaNames, cfg);

sessionTable = local_struct_to_table(sessionRows);
areaTable = local_struct_to_table(areaRows);

writetable(roleTable, fullfile(outDir, 'session_method_roles.csv'));
writetable(sessionTable, fullfile(outDir, 'session_summary.csv'));
writetable(areaTable, fullfile(outDir, 'session_area_summary.csv'));
save(fullfile(outDir, 'session_method_summary.mat'), ...
    'roleTable', 'sessionTable', 'areaTable', 'withinFile', 'crossFile', ...
    'planFile', 'finalFile', 'finalSessionIds', 'finalAreaNames');

fprintf('[summarize_ibl_method_sessions] sessions=%d  session-area rows=%d\n', ...
    height(sessionTable), height(areaTable));
fprintf('  wrote %s\n', fullfile(outDir, 'session_summary.csv'));
fprintf('  wrote %s\n', fullfile(outDir, 'session_area_summary.csv'));

function rows = local_collect_method_roles(withinFile, crossFile, planFile, finalFile)
rows = local_empty_role_rows();

if exist(withinFile, 'file')
    S = load(withinFile, 'val');
    if isfield(S, 'val') && isfield(S.val, 'records')
        recs = S.val.records;
        for i = 1:numel(recs)
            if isfield(recs(i), 'session_id')
                rows = local_add_role(rows, recs(i).session_id, ...
                    "step1_within_session", "target_full_session", ...
                    "within-session randomized subset validation");
            end
        end
    end
end

if exist(crossFile, 'file')
    S = load(crossFile, 'val');
    if isfield(S, 'val') && isfield(S.val, 'records')
        recs = S.val.records;
        for i = 1:numel(recs)
            if isfield(recs(i), 'target_session_id')
                rows = local_add_role(rows, recs(i).target_session_id, ...
                    "step2_cross_session", "target_full_session", ...
                    sprintf('cross-session candidate %d target', i));
            end
            if isfield(recs(i), 'sources')
                for s = 1:numel(recs(i).sources)
                    if isfield(recs(i).sources(s), 'session_id')
                        rows = local_add_role(rows, recs(i).sources(s).session_id, ...
                            "step2_cross_session", "source_partial_session", ...
                            sprintf('source for %s', string(recs(i).target_session_id)));
                    end
                end
            end
        end
    end
end

if exist(planFile, 'file')
    S = load(planFile, 'plan');
    if isfield(S, 'plan') && isfield(S.plan, 'selected_session_ids')
        ids = string(S.plan.selected_session_ids(:));
        for i = 1:numel(ids)
            rows = local_add_role(rows, ids(i), ...
                "step3_final_stitch", "planned_reliable_component", ...
                "selected by reliable-session component plan");
        end
    end
end

if exist(finalFile, 'file')
    S = load(finalFile, 'result');
    ids = local_result_kept_session_ids(S);
    for i = 1:numel(ids)
        rows = local_add_role(rows, ids(i), ...
            "step3_final_stitch", "actually_stitched", ...
            "kept in final stitched PDC result");
    end
end
end

function rows = local_empty_role_rows()
rows = struct('session_id', {}, 'step', {}, 'role', {}, 'detail', {});
end

function rows = local_add_role(rows, sessionId, step, role, detail)
if strlength(string(sessionId)) == 0
    return
end
row.session_id = string(sessionId);
row.step = string(step);
row.role = string(role);
row.detail = string(detail);
rows(end+1) = row; %#ok<AGROW>
end

function T = local_roles_to_table(rows)
if isempty(rows)
    T = table(strings(0,1), strings(0,1), strings(0,1), strings(0,1), ...
        'VariableNames', {'session_id', 'step', 'role', 'detail'});
    return
end
T = struct2table(rows);
end

function finalSessionIds = local_final_session_ids(planFile, finalFile)
finalSessionIds = strings(0,1);
if exist(finalFile, 'file')
    S = load(finalFile, 'result');
    finalSessionIds = local_result_kept_session_ids(S);
end
if isempty(finalSessionIds) && exist(planFile, 'file')
    S = load(planFile, 'plan');
    if isfield(S, 'plan') && isfield(S.plan, 'selected_session_ids')
        finalSessionIds = string(S.plan.selected_session_ids(:));
    end
end
finalSessionIds = unique(finalSessionIds(:), 'stable');
end

function ids = local_result_kept_session_ids(S)
ids = strings(0,1);
if ~isfield(S, 'result') || ~isfield(S.result, 'kept_sessions')
    return
end
kept = S.result.kept_sessions;
for i = 1:numel(kept)
    if isfield(kept(i), 'session_id')
        ids(end+1,1) = string(kept(i).session_id); %#ok<AGROW>
    end
end
end

function areaNames = local_final_area_names(finalFile)
areaNames = strings(0,1);
if exist(finalFile, 'file')
    S = load(finalFile, 'result');
    if isfield(S, 'result') && isfield(S.result, 'area_names')
        areaNames = string(S.result.area_names(:));
    end
end
end

function [sessionRows, areaRows] = local_build_summary_rows(sessionIds, roleTable, ...
    finalSessionIds, finalAreaNames, cfg)
sessionRows = struct([]);
areaRows = struct([]);

for s = 1:numel(sessionIds)
    sessionId = string(sessionIds(s));
    pcaFile = fullfile(cfg.area_pca_dir, sessionId + ".mat");
    [lab, subject, dateStr, numberStr] = local_parse_session_id(sessionId);
    rolesThis = roleTable(roleTable.session_id == sessionId, :);
    stepsInvolved = local_join_unique(rolesThis.step);
    rolesInvolved = local_join_unique(rolesThis.role);
    usedFinalSession = ismember(sessionId, finalSessionIds);

    if ~exist(pcaFile, 'file')
        row = local_session_row(sessionId, lab, subject, dateStr, numberStr, ...
            stepsInvolved, rolesInvolved, usedFinalSession, nan, nan, nan, ...
            strings(0,1), [], nan(0,1), "missing_pca_summary");
        sessionRows = local_append_struct(sessionRows, row);
        continue
    end

    S = load(pcaFile, 'summary');
    summary = S.summary;
    areaNames = string(summary.area_names(:));
    nNeurons = double(summary.n_neurons(:));
    durationS = local_get_scalar_field(summary, 'sp_dur', nan);
    pc1 = local_get_vector_field(summary, 'mean_pc1_explained', numel(areaNames), nan);
    passQc = local_get_logical_field(summary, 'pass_qc_area', numel(areaNames), true);
    qcReason = local_get_string_field(summary, 'qc_reason', numel(areaNames), "");

    row = local_session_row(sessionId, lab, subject, dateStr, numberStr, ...
        stepsInvolved, rolesInvolved, usedFinalSession, durationS, ...
        sum(nNeurons, 'omitnan'), numel(areaNames), areaNames, nNeurons, ...
        mean(pc1, 'omitnan'), "");
    row.n_areas_pass_qc = sum(passQc);
    row.area_neuron_counts = local_area_count_string(areaNames, nNeurons);
    row.areas_pass_qc = local_join_strings(areaNames(passQc));
    sessionRows = local_append_struct(sessionRows, row);

    for a = 1:numel(areaNames)
        areaUsedFinalNetwork = ismember(areaNames(a), finalAreaNames);
        sessionAreaUsedFinal = usedFinalSession && areaUsedFinalNetwork && passQc(a);
        reason = local_rejection_reason(usedFinalSession, areaUsedFinalNetwork, ...
            passQc(a), qcReason(a));

        arow.session_id = sessionId;
        arow.lab = lab;
        arow.subject = subject;
        arow.date = dateStr;
        arow.number = numberStr;
        arow.steps_involved = stepsInvolved;
        arow.roles_involved = rolesInvolved;
        arow.used_in_final_session = usedFinalSession;
        arow.area = areaNames(a);
        arow.n_neurons_area = nNeurons(a);
        arow.duration_s = durationS;
        arow.mean_pc1_explained = pc1(a);
        arow.pass_qc_area = passQc(a);
        arow.qc_reason = qcReason(a);
        arow.area_used_in_final_network = areaUsedFinalNetwork;
        arow.session_area_used_in_final = sessionAreaUsedFinal;
        arow.final_use_label = reason;
        areaRows = local_append_struct(areaRows, arow);
    end
end
end

function row = local_session_row(sessionId, lab, subject, dateStr, numberStr, ...
    stepsInvolved, rolesInvolved, usedFinalSession, durationS, nTotalNeurons, ...
    nAreas, areaNames, nNeurons, meanPc1, status)
row.session_id = sessionId;
row.lab = lab;
row.subject = subject;
row.date = dateStr;
row.number = numberStr;
row.steps_involved = stepsInvolved;
row.roles_involved = rolesInvolved;
row.used_in_final_session = usedFinalSession;
row.duration_s = durationS;
row.n_total_neurons = nTotalNeurons;
row.n_areas_total = nAreas;
row.n_areas_pass_qc = nan;
row.mean_pc1_explained = meanPc1;
row.areas_all = local_join_strings(areaNames);
row.areas_pass_qc = "";
row.area_neuron_counts = local_area_count_string(areaNames, nNeurons);
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

function [lab, subject, dateStr, numberStr] = local_parse_session_id(sessionId)
parts = split(string(sessionId), "__");
lab = "";
subject = "";
dateStr = "";
numberStr = "";
if numel(parts) >= 1, lab = parts(1); end
if numel(parts) >= 2, subject = parts(2); end
if numel(parts) >= 3, dateStr = parts(3); end
if numel(parts) >= 4, numberStr = parts(4); end
end

function x = local_get_scalar_field(S, name, defaultValue)
if isfield(S, name) && ~isempty(S.(name))
    x = double(S.(name));
else
    x = defaultValue;
end
end

function x = local_get_vector_field(S, name, n, defaultValue)
if isfield(S, name) && numel(S.(name)) == n
    value = S.(name);
    x = double(value(:));
else
    x = repmat(defaultValue, n, 1);
end
end

function x = local_get_logical_field(S, name, n, defaultValue)
if isfield(S, name) && numel(S.(name)) == n
    value = S.(name);
    x = logical(value(:));
else
    x = repmat(defaultValue, n, 1);
end
end

function x = local_get_string_field(S, name, n, defaultValue)
if isfield(S, name) && numel(S.(name)) == n
    value = S.(name);
    x = string(value(:));
else
    x = repmat(string(defaultValue), n, 1);
end
end

function s = local_join_unique(x)
x = unique(string(x(:)), 'stable');
s = local_join_strings(x);
end

function s = local_join_strings(x)
x = string(x(:));
x = x(strlength(x) > 0);
if isempty(x)
    s = "";
else
    s = strjoin(x.', "; ");
end
end

function s = local_area_count_string(areaNames, nNeurons)
areaNames = string(areaNames(:));
nNeurons = double(nNeurons(:));
if isempty(areaNames)
    s = "";
    return
end
parts = strings(numel(areaNames), 1);
for i = 1:numel(areaNames)
    if isnan(nNeurons(i))
        parts(i) = areaNames(i) + ":NA";
    else
        parts(i) = areaNames(i) + ":" + string(nNeurons(i));
    end
end
s = strjoin(parts.', "; ");
end

function reason = local_rejection_reason(usedFinalSession, areaUsedFinalNetwork, passQc, qcReason)
if usedFinalSession && areaUsedFinalNetwork && passQc
    reason = "used_final";
elseif ~passQc
    if strlength(string(qcReason)) > 0
        reason = "area_qc_failed: " + string(qcReason);
    else
        reason = "area_qc_failed";
    end
elseif ~usedFinalSession
    reason = "session_not_in_final_reliable_component";
elseif ~areaUsedFinalNetwork
    reason = "area_not_in_final_stitched_network";
else
    reason = "not_used_after_final_filters";
end
end
