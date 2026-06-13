function report = ibl_find_bad_session_spectra(cfg)
%IBL_FIND_BAD_SESSION_SPECTRA Identify cached session spectra with invalid explained variance.

files = dir(fullfile(cfg.session_spectra_dir, '*.mat'));

report = struct('file', {}, 'session_id', {}, 'n_areas', {}, ...
    'n_bad_scores', {}, 'bad_area_names', {}, 'has_nan_cross_spectrum', {}, ...
    'has_nan_auto_spectrum', {});

for k = 1:numel(files)
    pathk = fullfile(files(k).folder, files(k).name);
    tmp = load(pathk, 'summary');
    summary = tmp.summary;

    scores = local_scores(summary);
    badScores = ~isfinite(scores);
    badAreas = strings(0,1);
    if any(badScores)
        badAreas = string(summary.area_names(badScores));
    end

    hasNanCross = isfield(summary, 'cross_spectrum') && any(~isfinite(summary.cross_spectrum(:)));
    hasNanAuto = isfield(summary, 'auto_spectrum') && any(~isfinite(summary.auto_spectrum(:)));

    if any(badScores) || hasNanCross || hasNanAuto
        row.file = pathk;
        row.session_id = summary.session_id;
        row.n_areas = numel(summary.area_names);
        row.n_bad_scores = sum(badScores);
        row.bad_area_names = badAreas;
        row.has_nan_cross_spectrum = hasNanCross;
        row.has_nan_auto_spectrum = hasNanAuto;
        report(end+1) = row; %#ok<AGROW>
    end
end
end


function scores = local_scores(summary)
if isfield(summary, 'mean_pc1_explained')
    scores = summary.mean_pc1_explained(:);
else
    scores = nan(numel(summary.area_names), 1);
end

if ~isfield(summary, 'expl_var')
    return;
end

for i = 1:numel(scores)
    if isfinite(scores(i))
        continue;
    end
    vals = summary.expl_var{i};
    if isempty(vals)
        continue;
    end
    scores(i) = mean(vals(1,:), 'omitnan');
end
end
