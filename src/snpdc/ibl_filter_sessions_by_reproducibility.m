function filt = ibl_filter_sessions_by_reproducibility(cfg, varargin)
%IBL_FILTER_SESSIONS_BY_REPRODUCIBILITY Null-calibrated session state filter.
%
%   filt = ibl_filter_sessions_by_reproducibility(cfg)
%
% Scores each session by how well its observed pairwise CSD/coherence
% spectra match leave-one-session-out consensus spectra for the same area
% pairs. The default null preserves each session's number of pair
% observations but compares them to mismatched pair consensuses.
%
% The intended use is to remove whole sessions whose global coupling state
% is atypical, rather than discarding individual brain areas or pairs.

p = inputParser;
p.addParameter('Band', [1 80], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('Metric', 'coherence_abs', @(x)ischar(x)||isstring(x));
p.addParameter('MinSessionsPerPair', 5, @(x)isnumeric(x)&&isscalar(x)&&x>=3);
p.addParameter('MinPairsPerSession', 10, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('NumNull', 1000, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('Alpha', 0.05, @(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
p.addParameter('ThresholdMode', 'null', @(x)ischar(x)||isstring(x));
p.addParameter('RobustMadK', 2.5, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('Seed', 1, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('MakeFigure', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('OutputFile', "", @(x)isstring(x)||ischar(x));
p.parse(varargin{:});
opt = p.Results;
metric = lower(string(opt.Metric));
thresholdMode = lower(string(opt.ThresholdMode));

files = dir(fullfile(cfg.cross_spectra_dir, '*.mat'));
assert(~isempty(files), 'No cross-spectrum summaries found in %s.', cfg.cross_spectra_dir);

rows = local_load_pair_rows(cfg, files, opt, metric);
assert(~isempty(rows), 'No valid pair observations found in %s.', cfg.cross_spectra_dir);

[rows, pairInfo] = local_attach_consensus(rows, opt.MinSessionsPerPair);
assert(~isempty(rows), 'No pairs had at least MinSessionsPerPair=%d sessions.', opt.MinSessionsPerPair);

[sessionTable, pairObsTable] = local_score_sessions(rows, opt.MinPairsPerSession);
assert(~isempty(sessionTable), 'No sessions had at least MinPairsPerSession=%d scored pairs.', opt.MinPairsPerSession);

nullScores = local_null_session_scores(rows, sessionTable, pairInfo, opt);
[sessionTable, cutoff] = local_classify_sessions(sessionTable, nullScores, opt, thresholdMode);

filt.options = opt;
filt.metric = metric;
filt.threshold_mode = thresholdMode;
filt.score_cutoff = cutoff;
filt.session_table = sessionTable;
filt.pair_observation_table = pairObsTable;
filt.pair_table = local_pair_table(pairInfo);
filt.null_scores = nullScores;
filt.kept_session_ids = sessionTable.session_id(sessionTable.keep);
filt.excluded_session_ids = sessionTable.session_id(~sessionTable.keep);
filt.created_at = string(datetime('now'));

if strlength(string(opt.OutputFile)) > 0
    outFile = string(opt.OutputFile);
    outDir = fileparts(char(outFile));
    if ~isempty(outDir) && exist(outDir, 'dir') ~= 7
        mkdir(outDir);
    end
    save(outFile, 'filt', '-v7.3');
    fprintf('[ibl_filter_sessions_by_reproducibility] saved %s\n', char(outFile));
end

fprintf(['[ibl_filter_sessions_by_reproducibility] sessions=%d  kept=%d  excluded=%d  ' ...
    'cutoff=%.3f  median score=%.3f\n'], ...
    height(sessionTable), nnz(sessionTable.keep), nnz(~sessionTable.keep), ...
    cutoff, local_nanmedian(sessionTable.repro_score));

if opt.MakeFigure
    local_plot_filter(filt);
end
end


function rows = local_load_pair_rows(cfg, files, opt, metric)
rows = struct('session_id', {}, 'lab', {}, 'subject', {}, 'date', {}, ...
    'area_i', {}, 'area_j', {}, 'pair_key', {}, 'vector', {}, ...
    'band_mean', {});

fprintf('[ibl_filter_sessions_by_reproducibility] scanning %d session CSD files\n', numel(files));
for k = 1:numel(files)
    tmp = load(fullfile(files(k).folder, files(k).name), 'summary');
    summary = tmp.summary;
    if ~isfield(summary, 'cross_spectrum') || ~isfield(summary, 'freqs')
        continue;
    end

    freqs = summary.freqs(:);
    fMask = freqs >= opt.Band(1) & freqs <= min(opt.Band(2), freqs(end));
    if ~any(fMask)
        continue;
    end

    names = string(summary.area_names(:));
    parts = split(string(summary.session_id), "__");
    for i = 1:numel(names)
        for j = (i+1):numel(names)
            v = local_pair_vector(summary.cross_spectrum, i, j, fMask, metric);
            if isempty(v) || all(~isfinite(v))
                continue;
            end
            row.session_id = string(summary.session_id);
            row.lab = local_part(parts, 1);
            row.subject = local_part(parts, 2);
            row.date = local_part(parts, 3);
            row.area_i = names(i);
            row.area_j = names(j);
            row.pair_key = string(names(i) + "--" + names(j));
            row.vector = v(:).';
            row.band_mean = local_nanmean(v);
            rows(end+1) = row; %#ok<AGROW>
        end
    end

    if isfield(cfg, 'verbose') && cfg.verbose && (mod(k, 25) == 0 || k == numel(files))
        fprintf('  scanned %d/%d files, pair observations=%d\n', k, numel(files), numel(rows));
    end
end
end


function [rowsOut, pairInfo] = local_attach_consensus(rows, minSessionsPerPair)
pairKeys = unique(string({rows.pair_key}), 'stable');
rowsOut = struct('session_id', {}, 'lab', {}, 'subject', {}, 'date', {}, ...
    'area_i', {}, 'area_j', {}, 'pair_key', {}, 'vector', {}, ...
    'band_mean', {}, 'consensus_corr', {}, 'consensus_band_mean', {});
pairInfo = struct('pair_key', {}, 'area_i', {}, 'area_j', {}, ...
    'n_sessions', {}, 'session_ids', {}, 'vectors', {}, 'consensus', {});

for p = 1:numel(pairKeys)
    idx = find(string({rows.pair_key}) == pairKeys(p));
    sess = string({rows(idx).session_id});
    if numel(unique(sess)) < minSessionsPerPair
        continue;
    end

    vectors = vertcat(rows(idx).vector);
    consensus = local_nanmean_matrix(vectors);
    info.pair_key = pairKeys(p);
    info.area_i = rows(idx(1)).area_i;
    info.area_j = rows(idx(1)).area_j;
    info.n_sessions = numel(unique(sess));
    info.session_ids = sess(:);
    info.vectors = vectors;
    info.consensus = consensus;
    pairInfo(end+1) = info; %#ok<AGROW>

    for q = 1:numel(idx)
        thisRows = vectors;
        thisRows(q, :) = NaN;
        loo = local_nanmean_matrix(thisRows);
        r = local_vec_corr(rows(idx(q)).vector, loo);
        if ~isfinite(r)
            continue;
        end
        row = rows(idx(q));
        row.consensus_corr = r;
        row.consensus_band_mean = local_nanmean(loo);
        rowsOut(end+1) = row; %#ok<AGROW>
    end
end
end


function [sessionTable, pairObsTable] = local_score_sessions(rows, minPairsPerSession)
sessionIds = unique(string({rows.session_id}), 'stable');
sessionOut = struct('session_id', {}, 'lab', {}, 'subject', {}, 'date', {}, ...
    'n_pairs_scored', {}, 'repro_score', {}, 'mean_repro_score', {}, ...
    'iqr_low_repro_score', {}, 'iqr_high_repro_score', {}, ...
    'median_band_mean', {});

pairObsTable = local_pair_observation_table(rows);

for s = 1:numel(sessionIds)
    idx = find(string({rows.session_id}) == sessionIds(s));
    if numel(idx) < minPairsPerSession
        continue;
    end
    vals = [rows(idx).consensus_corr].';
    bandMeans = [rows(idx).band_mean].';
    out.session_id = sessionIds(s);
    out.lab = rows(idx(1)).lab;
    out.subject = rows(idx(1)).subject;
    out.date = rows(idx(1)).date;
    out.n_pairs_scored = numel(idx);
    out.repro_score = local_nanmedian(vals);
    out.mean_repro_score = local_nanmean(vals);
    out.iqr_low_repro_score = local_nanprctile(vals, 25);
    out.iqr_high_repro_score = local_nanprctile(vals, 75);
    out.median_band_mean = local_nanmedian(bandMeans);
    sessionOut(end+1) = out; %#ok<AGROW>
end

if isempty(sessionOut)
    sessionTable = table();
else
    sessionTable = struct2table(sessionOut);
end
end


function nullScores = local_null_session_scores(rows, sessionTable, pairInfo, opt)
rng(opt.Seed);
nSessions = height(sessionTable);
nullScores = nan(nSessions, opt.NumNull);
allPairKeys = string({pairInfo.pair_key});

fprintf('[ibl_filter_sessions_by_reproducibility] building null scores (%d shuffles)\n', opt.NumNull);
for s = 1:nSessions
    idx = find(string({rows.session_id}) == sessionTable.session_id(s));
    if isempty(idx)
        continue;
    end

    obsPairs = string({rows(idx).pair_key});
    obsVectors = {rows(idx).vector};
    for b = 1:opt.NumNull
        vals = nan(numel(idx), 1);
        for q = 1:numel(idx)
            candidate = find(allPairKeys ~= obsPairs(q));
            if isempty(candidate)
                continue;
            end
            pick = candidate(randi(numel(candidate)));
            vals(q) = local_vec_corr(obsVectors{q}, pairInfo(pick).consensus);
        end
        nullScores(s,b) = local_nanmedian(vals);
    end
    if mod(s, 25) == 0 || s == nSessions
        fprintf('  null session %d/%d\n', s, nSessions);
    end
end
end


function [T, cutoff] = local_classify_sessions(T, nullScores, opt, thresholdMode)
switch thresholdMode
    case "null"
        nullFlat = nullScores(isfinite(nullScores));
        cutoff = local_nanprctile(nullFlat, 100 * (1 - opt.Alpha));
        pRight = nan(height(T), 1);
        for i = 1:height(T)
            ns = nullScores(i, :);
            ns = ns(isfinite(ns));
            if isempty(ns)
                pRight(i) = NaN;
            else
                pRight(i) = (1 + nnz(ns >= T.repro_score(i))) / (1 + numel(ns));
            end
        end
        keep = pRight <= opt.Alpha;
        T.null_p_right = pRight;
        T.null_median = local_row_nanmedian(nullScores);
        T.null_p95 = local_row_nanprctile(nullScores, 95);
        T.keep_reason = strings(height(T),1);
        T.keep_reason(keep) = "above_null";
        T.keep_reason(~keep) = "not_above_null";
    case {"robust", "mad"}
        med = local_nanmedian(T.repro_score);
        madVal = local_mad(T.repro_score);
        cutoff = med - opt.RobustMadK * madVal;
        keep = T.repro_score >= cutoff;
        T.null_p_right = nan(height(T),1);
        T.null_median = local_row_nanmedian(nullScores);
        T.null_p95 = local_row_nanprctile(nullScores, 95);
        T.keep_reason = strings(height(T),1);
        T.keep_reason(keep) = "within_robust_range";
        T.keep_reason(~keep) = "low_tail_outlier";
    otherwise
        error('Unknown ThresholdMode: %s', thresholdMode);
end

T.keep = keep;
T = sortrows(T, {'keep', 'repro_score'}, {'descend', 'descend'});
end


function T = local_pair_observation_table(rows)
n = numel(rows);
session_id = strings(n,1);
lab = strings(n,1);
subject = strings(n,1);
date = strings(n,1);
pair_key = strings(n,1);
area_i = strings(n,1);
area_j = strings(n,1);
consensus_corr = nan(n,1);
band_mean = nan(n,1);
for i = 1:n
    session_id(i) = rows(i).session_id;
    lab(i) = rows(i).lab;
    subject(i) = rows(i).subject;
    date(i) = rows(i).date;
    pair_key(i) = rows(i).pair_key;
    area_i(i) = rows(i).area_i;
    area_j(i) = rows(i).area_j;
    consensus_corr(i) = rows(i).consensus_corr;
    band_mean(i) = rows(i).band_mean;
end
T = table(session_id, lab, subject, date, pair_key, area_i, area_j, ...
    consensus_corr, band_mean);
end


function T = local_pair_table(pairInfo)
if isempty(pairInfo)
    T = table();
    return;
end
n = numel(pairInfo);
pair_key = strings(n,1);
area_i = strings(n,1);
area_j = strings(n,1);
n_sessions = zeros(n,1);
for i = 1:n
    pair_key(i) = pairInfo(i).pair_key;
    area_i(i) = pairInfo(i).area_i;
    area_j(i) = pairInfo(i).area_j;
    n_sessions(i) = pairInfo(i).n_sessions;
end
T = table(pair_key, area_i, area_j, n_sessions);
T = sortrows(T, 'n_sessions', 'descend');
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


function y = local_nanmean_matrix(X)
if isempty(X)
    y = [];
    return;
end
y = nan(1, size(X,2));
for j = 1:size(X,2)
    y(j) = local_nanmean(X(:,j));
end
end


function y = local_row_nanmedian(X)
y = nan(size(X,1), 1);
for i = 1:size(X,1)
    y(i) = local_nanmedian(X(i,:));
end
end


function y = local_row_nanprctile(X, pct)
y = nan(size(X,1), 1);
for i = 1:size(X,1)
    y(i) = local_nanprctile(X(i,:), pct);
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


function y = local_nanprctile(x, pct)
x = x(isfinite(x));
if isempty(x)
    y = NaN;
else
    y = prctile(x, pct);
end
end


function y = local_mad(x)
x = x(isfinite(x));
if isempty(x)
    y = NaN;
else
    y = median(abs(x - median(x)));
end
end


function out = local_part(parts, idx)
if numel(parts) >= idx
    out = string(parts(idx));
else
    out = "";
end
end


function local_plot_filter(filt)
T = filt.session_table;
if isempty(T)
    return;
end
figure('Color', 'w', 'Position', [100 100 1120 430]);
tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
histogram(T.repro_score, 'BinEdges', -1:0.05:1, 'FaceColor', [0.2 0.45 0.75]);
hold on;
xline(filt.score_cutoff, 'r-', 'LineWidth', 1.5);
xlabel('Session reproducibility score');
ylabel('Session count');
title('Session-state filter');
grid on;

nexttile;
scatter(T.n_pairs_scored, T.repro_score, 42, double(T.keep), 'filled');
hold on;
yline(filt.score_cutoff, 'r-', 'LineWidth', 1.5);
xlabel('Scored pair observations');
ylabel('Session reproducibility score');
title('Coverage vs reproducibility');
colorbar;
grid on;
end
