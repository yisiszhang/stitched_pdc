function out = ibl_analyze_pair_coupling_variability(cfg, varargin)
%IBL_ANALYZE_PAIR_COUPLING_VARIABILITY Quantify CSD/coherence variability across sessions.
%
%   out = ibl_analyze_pair_coupling_variability(cfg)
%
% Reads cached cross-area CSD summaries and compares coupling spectra for
% the same area pair across sessions/animals. No Chronux spectra are
% recomputed.

p = inputParser;
p.addParameter('MinSessionsPerPair', 5, @(x)isnumeric(x)&&isscalar(x)&&x>=2);
p.addParameter('Band', [1 80], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('Metric', 'coherence_abs', @(x)ischar(x)||isstring(x));
p.addParameter('TopN', inf, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('MakeFigure', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('OutputFile', "", @(x)isstring(x)||ischar(x));
p.parse(varargin{:});
opt = p.Results;
metric = lower(string(opt.Metric));

files = dir(fullfile(cfg.cross_spectra_dir, '*.mat'));
assert(~isempty(files), 'No cross-spectrum summaries found in %s.', cfg.cross_spectra_dir);

sessionRows = struct('session_id', {}, 'lab', {}, 'subject', {}, 'date', {}, ...
    'area_i', {}, 'area_j', {}, 'pair_key', {}, 'n_freqs', {}, ...
    'band_mean', {}, 'band_max', {}, 'vector', {});

fprintf('[ibl_analyze_pair_coupling_variability] scanning %d session CSD files\n', numel(files));
for k = 1:numel(files)
    tmp = load(fullfile(files(k).folder, files(k).name), 'summary');
    summary = tmp.summary;
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
            row.n_freqs = numel(v);
            row.band_mean = local_nanmean(v);
            row.band_max = local_nanmax(v);
            row.vector = v(:).';
            sessionRows(end+1) = row; %#ok<AGROW>
        end
    end
    if isfield(cfg, 'verbose') && cfg.verbose && (mod(k, 25) == 0 || k == numel(files))
        fprintf('  scanned %d/%d files, pair observations=%d\n', k, numel(files), numel(sessionRows));
    end
end

sessionPairTable = local_session_pair_table(sessionRows);
pairTable = local_pair_variability_table(sessionRows, opt);

out.options = opt;
out.metric = metric;
out.session_pair_table = sessionPairTable;
out.pair_table = pairTable;
out.created_at = string(datetime('now'));

if strlength(string(opt.OutputFile)) > 0
    save(opt.OutputFile, 'out', '-v7.3');
end
if opt.MakeFigure
    local_plot_pair_variability(out);
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
        error('Unknown coupling variability metric: %s', metric);
end
v = real(v(:));
end


function T = local_session_pair_table(rows)
if isempty(rows)
    T = table();
    return;
end
n = numel(rows);
session_id = strings(n,1);
lab = strings(n,1);
subject = strings(n,1);
date = strings(n,1);
area_i = strings(n,1);
area_j = strings(n,1);
pair_key = strings(n,1);
n_freqs = zeros(n,1);
band_mean = nan(n,1);
band_max = nan(n,1);
for i = 1:n
    session_id(i) = rows(i).session_id;
    lab(i) = rows(i).lab;
    subject(i) = rows(i).subject;
    date(i) = rows(i).date;
    area_i(i) = rows(i).area_i;
    area_j(i) = rows(i).area_j;
    pair_key(i) = rows(i).pair_key;
    n_freqs(i) = rows(i).n_freqs;
    band_mean(i) = rows(i).band_mean;
    band_max(i) = rows(i).band_max;
end
T = table(session_id, lab, subject, date, area_i, area_j, pair_key, ...
    n_freqs, band_mean, band_max);
end


function T = local_pair_variability_table(rows, opt)
if isempty(rows)
    T = table();
    return;
end
pairKeys = unique(string({rows.pair_key}), 'stable');
outRows = struct('pair_key', {}, 'area_i', {}, 'area_j', {}, ...
    'n_sessions', {}, 'n_subjects', {}, 'n_labs', {}, ...
    'median_session_corr', {}, 'mean_session_corr', {}, ...
    'iqr_low_session_corr', {}, 'iqr_high_session_corr', {}, ...
    'frac_corr_gt_0p5', {}, 'band_mean_median', {}, 'band_mean_cv', {}, ...
    'within_lab_median_corr', {}, 'across_lab_median_corr', {});

for p = 1:numel(pairKeys)
    idx = find(string({rows.pair_key}) == pairKeys(p));
    if numel(idx) < opt.MinSessionsPerPair
        continue;
    end
    corrs = [];
    withinLabCorrs = [];
    acrossLabCorrs = [];
    for a = 1:numel(idx)
        for b = (a+1):numel(idx)
            ia = idx(a);
            ib = idx(b);
            r = local_vec_corr(rows(ia).vector, rows(ib).vector);
            corrs(end+1,1) = r; %#ok<AGROW>
            if rows(ia).lab == rows(ib).lab
                withinLabCorrs(end+1,1) = r; %#ok<AGROW>
            else
                acrossLabCorrs(end+1,1) = r; %#ok<AGROW>
            end
        end
    end
    bandMeans = [rows(idx).band_mean].';
    row.pair_key = pairKeys(p);
    row.area_i = rows(idx(1)).area_i;
    row.area_j = rows(idx(1)).area_j;
    row.n_sessions = numel(idx);
    row.n_subjects = numel(unique(string({rows(idx).subject})));
    row.n_labs = numel(unique(string({rows(idx).lab})));
    row.median_session_corr = local_nanmedian(corrs);
    row.mean_session_corr = local_nanmean(corrs);
    row.iqr_low_session_corr = local_nanprctile(corrs, 25);
    row.iqr_high_session_corr = local_nanprctile(corrs, 75);
    row.frac_corr_gt_0p5 = local_nanmean(double(corrs > 0.5));
    row.band_mean_median = local_nanmedian(bandMeans);
    row.band_mean_cv = local_nanstd(bandMeans) / max(abs(local_nanmean(bandMeans)), eps);
    row.within_lab_median_corr = local_nanmedian(withinLabCorrs);
    row.across_lab_median_corr = local_nanmedian(acrossLabCorrs);
    outRows(end+1) = row; %#ok<AGROW>
end

if isempty(outRows)
    T = table();
    return;
end
T = struct2table(outRows);
T = sortrows(T, {'n_sessions', 'median_session_corr'}, {'descend', 'descend'});
if isfinite(opt.TopN) && height(T) > opt.TopN
    T = T(1:opt.TopN, :);
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


function y = local_nanmax(x)
x = x(isfinite(x));
if isempty(x)
    y = NaN;
else
    y = max(x);
end
end


function y = local_nanstd(x)
x = x(isfinite(x));
if numel(x) < 2
    y = NaN;
else
    y = std(x);
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


function out = local_part(parts, idx)
if numel(parts) >= idx
    out = string(parts(idx));
else
    out = "";
end
end


function local_plot_pair_variability(out)
T = out.pair_table;
if isempty(T)
    return;
end
figure('Color', 'w', 'Position', [100 100 1120 430]);
tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
histogram(T.median_session_corr, 'BinEdges', -1:0.1:1);
xlabel('Median session-to-session coupling-spectrum corr');
ylabel('Area-pair count');
title(sprintf('Pair coupling reproducibility (%s)', out.metric), 'Interpreter', 'none');
grid on;

nexttile;
scatter(T.n_sessions, T.median_session_corr, 45, T.n_labs, 'filled');
xlabel('Number of sessions observing pair');
ylabel('Median coupling-spectrum corr');
title('Reproducibility vs support');
colorbar;
grid on;
end
