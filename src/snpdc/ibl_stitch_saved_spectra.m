function result = ibl_stitch_saved_spectra(cfg, scan)
%IBL_STITCH_SAVED_SPECTRA Stitch per-session area spectra and compute PDC.

if nargin < 2 || isempty(scan)
    if exist(cfg.pca_qc_file, 'file') == 2
        tmp = load(cfg.pca_qc_file, 'qc');
        scan = tmp.qc;
    else
        tmp = load(cfg.scan_file, 'scan');
        scan = tmp.scan;
    end
end

if isfield(scan, 'largest_component') && isfield(scan, 'qualifying_sessions')
    globalAreas = scan.largest_component(:);
else
    globalAreas = scan.largest_component(:);
end

if exist(cfg.cross_spectra_dir, 'dir') == 7
    files = dir(fullfile(cfg.cross_spectra_dir, '*.mat'));
else
    files = dir(fullfile(cfg.session_spectra_dir, '*.mat'));
end
assert(~isempty(files), 'No session spectra found in %s', cfg.session_spectra_dir);

if isfield(cfg, 'verbose') && cfg.verbose
    fprintf('[ibl_stitch_saved_spectra] Found %d saved session spectra\n', numel(files));
end

sessionKeep = local_session_keep_set(cfg);

Sblocks = {};
recset = {};
kept = struct('session_id', {}, 'areas', {});
f = [];

for k = 1:numel(files)
    tmp = load(fullfile(files(k).folder, files(k).name), 'summary');
    summary = tmp.summary;
    if isfield(cfg, 'verbose') && cfg.verbose
        fprintf('[ibl_stitch_saved_spectra] %4d/%4d  loading %s\n', ...
            k, numel(files), summary.session_id);
    end
    if ~local_keep_session(summary.session_id, sessionKeep)
        if isfield(cfg, 'verbose') && cfg.verbose
            fprintf('  skip: excluded by session filter\n');
        end
        continue;
    end

    areaScores = local_area_scores(summary);
    keepArea = areaScores >= cfg.min_pc1_explained;
    areaNames = string(summary.area_names(keepArea));
    if numel(areaNames) < 2
        if isfield(cfg, 'verbose') && cfg.verbose
            fprintf('  skip: only %d areas after explained-variance threshold\n', numel(areaNames));
        end
        continue;
    end

    [tf, idxGlobal] = ismember(areaNames, globalAreas);
    areaNames = areaNames(tf);
    idxGlobal = idxGlobal(tf);
    if numel(idxGlobal) < 2
        if isfield(cfg, 'verbose') && cfg.verbose
            fprintf('  skip: only %d areas overlap with global component\n', numel(idxGlobal));
        end
        continue;
    end

    block = summary.cross_spectrum(keepArea, keepArea, :);
    block = block(tf, tf, :);

    if isempty(f)
        f = summary.freqs(:);
    else
        assert(numel(f) == numel(summary.freqs) && max(abs(f - summary.freqs(:))) < 1e-10, ...
            'Frequency mismatch across session summaries.');
    end

    Sblocks{end+1} = block; %#ok<AGROW>
    recset{end+1} = idxGlobal(:)'; %#ok<AGROW>
    kept(end+1).session_id = summary.session_id; %#ok<AGROW>
    kept(end).areas = areaNames;
    if isfield(cfg, 'verbose') && cfg.verbose
        fprintf('  keep: %d areas\n', numel(areaNames));
    end
end

assert(~isempty(Sblocks), 'No saved session spectra passed the stitching filters.');

stitchParams = cfg.stitch;
stitchParams.verbose = isfield(cfg, 'verbose') && cfg.verbose;
[P, S, Sreg, f, meta] = stitch_spectra_blocks(Sblocks, recset, f, stitchParams);
if isfield(cfg, 'verbose') && cfg.verbose
    fprintf('[ibl_stitch_saved_spectra] running Wilson-factorization PDC\n');
end
PDC = nonparam_pdc_H(Sreg, f, 'maxiter', cfg.chronux.maxiter, 'tol', cfg.chronux.tol);
if isfield(cfg, 'verbose') && cfg.verbose
    fprintf('[ibl_stitch_saved_spectra] done\n');
end

result.area_names = globalAreas;
result.freqs = f;
result.S = S;
result.Sreg = Sreg;
result.P = P;
result.PDC = PDC;
result.recset = recset;
result.kept_sessions = kept;
result.meta = meta;

save(cfg.stitch_file, 'result', '-v7.3');
end


function sessionKeep = local_session_keep_set(cfg)
sessionKeep = [];
if isfield(cfg, 'session_include_ids') && ~isempty(cfg.session_include_ids)
    sessionKeep.include = string(cfg.session_include_ids(:));
else
    sessionKeep.include = strings(0,1);
end
if isfield(cfg, 'session_exclude_ids') && ~isempty(cfg.session_exclude_ids)
    sessionKeep.exclude = string(cfg.session_exclude_ids(:));
else
    sessionKeep.exclude = strings(0,1);
end

if isfield(cfg, 'session_filter_file') && strlength(string(cfg.session_filter_file)) > 0
    filterFile = string(cfg.session_filter_file);
    assert(exist(filterFile, 'file') == 2, 'session_filter_file not found: %s', char(filterFile));
    tmp = load(filterFile, 'filt');
    assert(isfield(tmp, 'filt') && isfield(tmp.filt, 'kept_session_ids'), ...
        'session_filter_file must contain filt.kept_session_ids.');
    filterInclude = string(tmp.filt.kept_session_ids(:));
    if isempty(sessionKeep.include)
        sessionKeep.include = filterInclude;
    else
        sessionKeep.include = intersect(sessionKeep.include, filterInclude, 'stable');
    end
end
end


function tf = local_keep_session(sessionId, sessionKeep)
if isempty(sessionKeep)
    tf = true;
    return;
end
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


function scores = local_area_scores(summary)
if isfield(summary, 'mean_pc1_explained')
    scores = summary.mean_pc1_explained(:);
else
    scores = nan(numel(summary.area_names), 1);
end

bad = ~isfinite(scores);
if ~any(bad)
    return;
end

if ~isfield(summary, 'expl_var') || ~isfield(summary, 'freqs')
    scores(bad) = -inf;
    return;
end

freqs = summary.freqs(:);
if isfield(summary, 'pc_band') && numel(summary.pc_band) == 2
    band = summary.pc_band;
else
    band = [1 100];
end
mask = freqs >= band(1) & freqs <= min(band(2), freqs(end));
if ~any(mask)
    mask = true(size(freqs));
end

for i = 1:numel(scores)
    if ~bad(i)
        continue;
    end
    if i > numel(summary.expl_var) || isempty(summary.expl_var{i})
        scores(i) = -inf;
        continue;
    end
    vals = summary.expl_var{i};
    if isempty(vals) || size(vals,1) < 1
        scores(i) = -inf;
        continue;
    end
    scores(i) = mean(vals(1, mask), 'omitnan');
    if ~isfinite(scores(i))
        scores(i) = -inf;
    end
end
end
