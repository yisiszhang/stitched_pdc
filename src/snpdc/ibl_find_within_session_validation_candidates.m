function candidates = ibl_find_within_session_validation_candidates(qc, cfg, varargin)
%IBL_FIND_WITHIN_SESSION_VALIDATION_CANDIDATES Rank sessions for split-half validation.
%
%   candidates = ibl_find_within_session_validation_candidates(qc, cfg)
%
% This is a fast screening step. It uses cached area-PCA summaries and does
% not recompute Chronux spectra. A good candidate has long spontaneous
% duration and many PCA-QC-passed areas inside the current largest component.

p = inputParser;
p.addParameter('MinAreas', 6, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('MaxAreas', inf, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('MinDuration', 600, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('MinMeanPc1', [], @(x)isnumeric(x)||isempty(x));
p.addParameter('TopN', 30, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('RequirePcaCache', true, @(x)islogical(x)&&isscalar(x));
p.parse(varargin{:});
opt = p.Results;

if nargin < 1 || isempty(qc)
    tmp = load(cfg.pca_qc_file, 'qc');
    qc = tmp.qc;
end
if isempty(opt.MinMeanPc1)
    opt.MinMeanPc1 = cfg.min_pc1_explained;
end

sessions = qc.qualifying_sessions;
rows = struct('session_id', {}, 'lab', {}, 'subject', {}, 'date', {}, ...
    'duration_s', {}, 'n_pass_areas', {}, 'n_selected_areas', {}, ...
    'mean_pc1', {}, 'median_pc1', {}, 'min_pc1', {}, ...
    'area_names', {}, 'pca_file', {});

for k = 1:numel(sessions)
    sessionId = string(sessions(k).session_id);
    pcaFile = fullfile(cfg.area_pca_dir, char(sessionId + ".mat"));
    if exist(pcaFile, 'file') ~= 2
        if opt.RequirePcaCache
            continue;
        end
        spDur = sessions(k).sp_dur;
        areaNames = string(sessions(k).area_names(:));
        scores = nan(size(areaNames));
    else
        tmp = load(pcaFile, 'summary');
        summary = tmp.summary;
        spDur = summary.sp_dur;
        scoresAll = summary.mean_pc1_explained(:);
        keep = summary.pass_qc_area(:) & ...
            ismember(string(summary.area_names(:)), qc.largest_component(:)) & ...
            isfinite(scoresAll) & scoresAll >= opt.MinMeanPc1;
        areaNames = string(summary.area_names(keep));
        scores = scoresAll(keep);
    end

    if spDur < opt.MinDuration || numel(areaNames) < opt.MinAreas
        continue;
    end

    [scores, ord] = sort(scores, 'descend');
    areaNames = areaNames(ord);
    nPassAreas = numel(areaNames);
    if isfinite(opt.MaxAreas)
        take = 1:min(numel(areaNames), opt.MaxAreas);
        areaNames = areaNames(take);
        scores = scores(take);
    end

    parts = split(sessionId, "__");
    row.session_id = sessionId;
    row.lab = local_part(parts, 1);
    row.subject = local_part(parts, 2);
    row.date = local_part(parts, 3);
    row.duration_s = spDur;
    row.n_pass_areas = nPassAreas;
    row.n_selected_areas = numel(areaNames);
    row.mean_pc1 = mean(scores, 'omitnan');
    row.median_pc1 = median(scores, 'omitnan');
    row.min_pc1 = min(scores, [], 'omitnan');
    row.area_names = strjoin(areaNames(:).', ", ");
    row.pca_file = string(pcaFile);
    rows(end+1) = row; %#ok<AGROW>
end

if isempty(rows)
    candidates = table();
    return;
end

candidates = struct2table(rows);
candidates = sortrows(candidates, ...
    {'n_pass_areas', 'duration_s', 'mean_pc1'}, ...
    {'descend', 'descend', 'descend'});
if isfinite(opt.TopN) && height(candidates) > opt.TopN
    candidates = candidates(1:opt.TopN, :);
end
end


function out = local_part(parts, idx)
if numel(parts) >= idx
    out = string(parts(idx));
else
    out = "";
end
end
