function diagOut = ibl_diagnose_pair_fingerprint_vs_csd(valCross, fp, varargin)
%IBL_DIAGNOSE_PAIR_FINGERPRINT_VS_CSD Relate area fingerprints to pair CSD.
%
%   diagOut = ibl_diagnose_pair_fingerprint_vs_csd(valCross, fp)
%
% Uses a cross-session validation result and spectral-fingerprint diagnostic
% for the same candidate. For each observed source-supported area pair, it
% compares source-stitched vs target CSD similarity to the source-target
% fingerprint similarity of the two endpoint areas.

p = inputParser;
p.addParameter('RecordIndex', 1, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('FingerprintMode', 'max', @(x)ischar(x)||isstring(x));
p.addParameter('MakeFigure', true, @(x)islogical(x)&&isscalar(x));
p.addParameter('OutputFile', "", @(x)isstring(x)||ischar(x));
p.parse(varargin{:});
opt = p.Results;
mode = lower(string(opt.FingerprintMode));

if isfield(valCross, 'records')
    rec = valCross.records(opt.RecordIndex);
else
    rec = valCross(opt.RecordIndex);
end

areaNames = string(rec.covered_area_names(:));
Sstitch = rec.S_stitched_sources;
Starget = rec.S_target_full_control_metric;
observedPair = ~rec.missing_mask;
completedPair = rec.missing_mask;
off = ~eye(numel(areaNames));

rows = struct('area_i', {}, 'area_j', {}, 'pair_type', {}, ...
    'csd_abs_corr', {}, 'csd_rel_rmse', {}, ...
    'infoflow_stitch_ij', {}, 'infoflow_target_ij', {}, ...
    'infoflow_stitch_ji', {}, 'infoflow_target_ji', {}, ...
    'fingerprint_i', {}, 'fingerprint_j', {}, 'pair_fingerprint', {}, ...
    'explained_var_i', {}, 'explained_var_j', {}, ...
    'eig_spectrum_i', {}, 'eig_spectrum_j', {}, ...
    'power_i', {}, 'power_j', {});

for i = 1:numel(areaNames)
    for j = (i+1):numel(areaNames)
        if ~off(i,j)
            continue;
        end
        row.area_i = areaNames(i);
        row.area_j = areaNames(j);
        if observedPair(i,j)
            row.pair_type = "observed";
        elseif completedPair(i,j)
            row.pair_type = "completed";
        else
            row.pair_type = "other";
        end
        row.csd_abs_corr = local_vec_corr(abs(squeeze(Sstitch(i,j,:))), ...
            abs(squeeze(Starget(i,j,:))));
        row.csd_rel_rmse = norm(squeeze(Sstitch(i,j,:)) - squeeze(Starget(i,j,:))) / ...
            max(norm(squeeze(Starget(i,j,:))), eps);
        if isfield(rec, 'infoflow_stitched_sources')
            row.infoflow_stitch_ij = rec.infoflow_stitched_sources(i,j);
            row.infoflow_target_ij = rec.infoflow_target_full_control(i,j);
            row.infoflow_stitch_ji = rec.infoflow_stitched_sources(j,i);
            row.infoflow_target_ji = rec.infoflow_target_full_control(j,i);
        else
            row.infoflow_stitch_ij = NaN;
            row.infoflow_target_ij = NaN;
            row.infoflow_stitch_ji = NaN;
            row.infoflow_target_ji = NaN;
        end

        statI = local_area_fp_stats(fp.area_table, areaNames(i), mode);
        statJ = local_area_fp_stats(fp.area_table, areaNames(j), mode);
        row.fingerprint_i = statI.feature_corr;
        row.fingerprint_j = statJ.feature_corr;
        row.pair_fingerprint = local_nanmean([row.fingerprint_i row.fingerprint_j]);
        row.explained_var_i = statI.explained_var_corr;
        row.explained_var_j = statJ.explained_var_corr;
        row.eig_spectrum_i = statI.eig_spectrum_corr;
        row.eig_spectrum_j = statJ.eig_spectrum_corr;
        row.power_i = statI.power_corr;
        row.power_j = statJ.power_corr;
        rows(end+1) = row; %#ok<AGROW>
    end
end

pairTable = struct2table(rows);
observedMask = pairTable.pair_type == "observed";
completedMask = pairTable.pair_type == "completed";

summary = struct();
summary.n_pairs = height(pairTable);
summary.n_observed_pairs = nnz(observedMask);
summary.n_completed_pairs = nnz(completedMask);
summary.observed_pair_fingerprint_vs_csd_corr = local_vec_corr( ...
    pairTable.pair_fingerprint(observedMask), pairTable.csd_abs_corr(observedMask));
summary.completed_pair_fingerprint_vs_csd_corr = local_vec_corr( ...
    pairTable.pair_fingerprint(completedMask), pairTable.csd_abs_corr(completedMask));
summary.observed_median_csd_corr_high_fp = local_median_by_threshold( ...
    pairTable.csd_abs_corr(observedMask), pairTable.pair_fingerprint(observedMask), 0.8, true);
summary.observed_median_csd_corr_low_fp = local_median_by_threshold( ...
    pairTable.csd_abs_corr(observedMask), pairTable.pair_fingerprint(observedMask), 0.8, false);

diagOut.record_index = opt.RecordIndex;
diagOut.target_session_id = string(rec.target_session_id);
diagOut.fingerprint_mode = mode;
diagOut.pair_table = pairTable;
diagOut.summary = summary;
diagOut.created_at = string(datetime('now'));

if strlength(string(opt.OutputFile)) > 0
    save(opt.OutputFile, 'diagOut', '-v7.3');
end
if opt.MakeFigure
    local_plot_pair_diag(diagOut);
end
end


function stat = local_area_fp_stats(areaTable, area, mode)
mask = areaTable.area == area;
if ~any(mask)
    stat.feature_corr = NaN;
    stat.explained_var_corr = NaN;
    stat.eig_spectrum_corr = NaN;
    stat.power_corr = NaN;
    return;
end
switch mode
    case "median"
        stat.feature_corr = local_nanmedian(areaTable.feature_corr(mask));
        stat.explained_var_corr = local_nanmedian(areaTable.explained_var_corr(mask));
        stat.eig_spectrum_corr = local_nanmedian(areaTable.eig_spectrum_corr(mask));
        stat.power_corr = local_nanmedian(areaTable.power_corr(mask));
    otherwise
        stat.feature_corr = local_nanmax(areaTable.feature_corr(mask));
        stat.explained_var_corr = local_nanmax(areaTable.explained_var_corr(mask));
        stat.eig_spectrum_corr = local_nanmax(areaTable.eig_spectrum_corr(mask));
        stat.power_corr = local_nanmax(areaTable.power_corr(mask));
end
end


function y = local_median_by_threshold(values, scores, thresh, high)
if high
    mask = scores >= thresh;
else
    mask = scores < thresh;
end
values = values(mask);
y = local_nanmedian(values);
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


function local_plot_pair_diag(diagOut)
T = diagOut.pair_table;
figure('Color', 'w', 'Position', [100 100 1120 430]);
tiledlayout(1,2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
hold on;
obs = T.pair_type == "observed";
comp = T.pair_type == "completed";
scatter(T.pair_fingerprint(obs), T.csd_abs_corr(obs), 55, 'filled', 'MarkerFaceAlpha', 0.75);
scatter(T.pair_fingerprint(comp), T.csd_abs_corr(comp), 55, 'filled', 'MarkerFaceAlpha', 0.75);
xlabel('Mean endpoint spectral-fingerprint similarity');
ylabel('Pair CSD magnitude correlation');
legend({'observed', 'completed'}, 'Location', 'best');
title(sprintf('Target %s', diagOut.target_session_id), 'Interpreter', 'none');
grid on;

nexttile;
hold on;
scatter(T.pair_fingerprint(obs), T.csd_rel_rmse(obs), 55, 'filled', 'MarkerFaceAlpha', 0.75);
scatter(T.pair_fingerprint(comp), T.csd_rel_rmse(comp), 55, 'filled', 'MarkerFaceAlpha', 0.75);
xlabel('Mean endpoint spectral-fingerprint similarity');
ylabel('Pair CSD relative RMSE');
legend({'observed', 'completed'}, 'Location', 'best');
title('Scale/error vs marginal similarity');
grid on;
end
