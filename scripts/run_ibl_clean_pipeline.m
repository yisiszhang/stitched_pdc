%RUN_IBL_CLEAN_PIPELINE Fresh end-to-end IBL stitched-PDC pipeline.
%
% This script starts from raw cached IBL files and writes every derived file
% into a new run folder under ibl_output/clean_runs. It does not read old
% area_pca, cross_spectra, pca_qc, stitch, or validation outputs.
%
% Step outputs are intentionally explicit:
%   00_scan/session_scan.mat                  -> raw-session scan
%   01_area_pca/*.mat                         -> per-session/area spectral PCA
%   02_pca_qc/pca_qc.mat                      -> PCA-QC graph/component
%   02_within_session_validation/*.mat,*.csv   -> split-half subset validation
%   03_cross_spectra/*.mat                    -> CSD projected with matching PCs
%   04_session_filter/session_repro_filter.mat-> optional session-state filter
%   05_reliable_plan/reliable_plan.mat        -> selected reliable component
%   06_stitch/stitched_pdc.mat                -> stitched S, PDC, metadata
%   07_infoflow_allen/*.mat,*.csv             -> info-flow and Allen comparison

clearvars;
setup;

% -------------------------------------------------------------------------
% Transparent configuration. Edit only this block for a new analysis.
% -------------------------------------------------------------------------
runTag = "clean_maxneurons100";
maxNeuronsPerArea = 100;
targetNFreqs = 128;
pcBand = [1 100];
analysisBand = [1 80];

minSpDur = 300;
minMeanFr = 0.2;
minNeuronsPerArea = 5;
minPc1Explained = 0.10;
minSessionsPerArea = 3;

runWithinSessionValidation = true;
withinValidationNumRepeats = 20;
withinValidationMaxSessions = 2;
withinValidationMinAreas = 6;
withinValidationMaxAreas = 12;
withinValidationMinDuration = 600;
withinValidationNumObservations = 2;
withinValidationOverlapFraction = 0.35;
withinValidationSeed = 1;
withinValidationParallel = true;

applySessionReproFilter = true;
reproMinSessionsPerPair = 3;
reproMinPairsPerSession = 5;
reproThresholdMode = "robust";
reproRobustMadK = 2.5;
reproNumNull = 1000;

deriveCsdCorrThreshold = true;
csdCorrThresholdFallback = 0.30;
thresholdOutcome = "info_observed_corr"; % lower 95% CI of this vs observed CSD corr crosses 0
crossValidationTopN = 20;
crossValidationMaxCandidates = 10;
crossValidationMinTargetAreas = 6;
crossValidationMaxTargetAreas = 14;
crossValidationMinSourceAreas = 3;
crossValidationMinDuration = 300;
crossValidationMinCoverage = 0.5;
crossValidationMinAnchorAreas = 2;
crossValidationMinAnchorFraction = 0.0;
crossValidationMinObservedPairFraction = 0.0;
crossValidationMaxCompletedPairFraction = 1.0;
crossValidationRequireConnectedSources = true;
crossValidationMinSourceOverlapAreas = 1;

planMinOverlapAreas = 2;
planMinOverlapPairs = 1;
planOverlapComponentSelection = "largest_sessions";
planReliableComponentSelection = "largest_areas";

stitchNormalize = "coherence";  % "none" or "coherence"
stitchRegularizer = "eigfloor"; % "eigfloor" or "glasso"
stitchLambda = 0;
stitchParallel = false;

allenMatrixCsv = fullfile(pwd, 'ibl_output', 'allen_tracing', ...
    'allen_tracing_matrix__normalized_projection_volume__max.csv');
allowExistingRunFolder = false;

% -------------------------------------------------------------------------
% Run folder and config paths.
% -------------------------------------------------------------------------
runStamp = string(datestr(now, 'yyyymmdd_HHMMSS'));
runName = runTag + "_" + runStamp;
baseCfg = ibl_default_config();
runRoot = fullfile(baseCfg.output_root, 'clean_runs', char(runName));
if exist(runRoot, 'dir') == 7 && ~allowExistingRunFolder
    error('Run folder already exists: %s', runRoot);
end

cfg = ibl_default_config(baseCfg.cache_root, runRoot);
cfg = local_apply_clean_paths(cfg, runRoot);
cfg.max_neurons_per_area = maxNeuronsPerArea;
cfg.target_n_freqs = targetNFreqs;
cfg.pc_band = pcBand;
cfg.min_sp_dur = minSpDur;
cfg.min_mean_fr = minMeanFr;
cfg.min_neurons_per_area = minNeuronsPerArea;
cfg.min_pc1_explained = minPc1Explained;
cfg.min_sessions_per_area = minSessionsPerArea;
cfg.verbose = true;
cfg.progress_every = 10;
cfg.stitch.normalize = char(stitchNormalize);
cfg.stitch.regularizer = char(stitchRegularizer);
cfg.stitch.lambda = stitchLambda;
cfg.stitch.parallel = stitchParallel;
cfg.threshold_validation = struct( ...
    'derive_csd_corr_threshold', deriveCsdCorrThreshold, ...
    'fallback_threshold', csdCorrThresholdFallback, ...
    'outcome', thresholdOutcome, ...
    'top_n', crossValidationTopN, ...
    'max_candidates', crossValidationMaxCandidates, ...
    'min_target_areas', crossValidationMinTargetAreas, ...
    'max_target_areas', crossValidationMaxTargetAreas, ...
    'min_source_areas', crossValidationMinSourceAreas, ...
    'min_duration', crossValidationMinDuration, ...
    'min_coverage', crossValidationMinCoverage, ...
    'min_anchor_areas', crossValidationMinAnchorAreas, ...
    'min_anchor_fraction', crossValidationMinAnchorFraction, ...
    'min_observed_pair_fraction', crossValidationMinObservedPairFraction, ...
    'max_completed_pair_fraction', crossValidationMaxCompletedPairFraction, ...
    'require_connected_sources', crossValidationRequireConnectedSources, ...
    'min_source_overlap_areas', crossValidationMinSourceOverlapAreas);
cfg.within_session_validation = struct( ...
    'enabled', runWithinSessionValidation, ...
    'num_repeats', withinValidationNumRepeats, ...
    'max_sessions', withinValidationMaxSessions, ...
    'min_areas', withinValidationMinAreas, ...
    'max_areas', withinValidationMaxAreas, ...
    'min_duration', withinValidationMinDuration, ...
    'num_observations', withinValidationNumObservations, ...
    'overlap_fraction', withinValidationOverlapFraction, ...
    'seed', withinValidationSeed, ...
    'parallel', withinValidationParallel);

local_make_clean_dirs(cfg);
local_write_config_summary(cfg, runRoot, runName, analysisBand, ...
    applySessionReproFilter, reproMinSessionsPerPair, reproMinPairsPerSession, ...
    reproThresholdMode, reproRobustMadK, reproNumNull, deriveCsdCorrThreshold, ...
    csdCorrThresholdFallback, thresholdOutcome, ...
    planMinOverlapAreas, planMinOverlapPairs, planOverlapComponentSelection, ...
    planReliableComponentSelection, allenMatrixCsv);

fprintf('[clean pipeline] run folder: %s\n', runRoot);

%%
% -------------------------------------------------------------------------
% Step 0: scan raw sessions.
% Saves cfg.scan_file. Used by Step 1.
% -------------------------------------------------------------------------
fprintf('\n[clean pipeline] Step 0/7: scan raw sessions\n');
scan = ibl_scan_sessions(cfg);
local_write_scan_summary(scan, fullfile(runRoot, '00_scan'));

%%
% -------------------------------------------------------------------------
% Step 1: per-session/area spectral PCA.
% Saves one .mat per session in cfg.area_pca_dir. Used by Step 2 and Step 3.
% -------------------------------------------------------------------------
fprintf('\n[clean pipeline] Step 1/7: spectral PCA per session/area\n');
ibl_compute_area_pca(scan, cfg);
local_write_area_pca_summary(cfg, fullfile(runRoot, '01_area_pca'));

%%
% -------------------------------------------------------------------------
% Step 2: PCA-QC graph/component.
% Saves cfg.pca_qc_file. Used by Step 2b, Step 3, and final stitch universe.
% -------------------------------------------------------------------------
fprintf('\n[clean pipeline] Step 2/7: PCA QC graph\n');
qc = ibl_build_pca_qc_graph(cfg, scan);
local_write_qc_summary(qc, fullfile(runRoot, '02_pca_qc'));

%%
% -------------------------------------------------------------------------
% Step 2b: within-session split-half/randomized-subset validation.
% Saves validation records and summaries. This is a method check only; later
% CSD/stitching steps do not read these files. 
% -------------------------------------------------------------------------
fprintf('\n[clean pipeline] Step 2b/7: within-session split-half subset validation\n');
withinDir = fullfile(runRoot, '02_within_session_validation');
if runWithinSessionValidation
    withinCandidates = ibl_find_within_session_validation_candidates(qc, cfg, ...
        'MinAreas', withinValidationMinAreas, ...
        'MaxAreas', withinValidationMaxAreas, ...
        'MinDuration', withinValidationMinDuration, ...
        'TopN', 30);
    writetable(withinCandidates, fullfile(withinDir, 'within_session_candidates.csv'));
    if isempty(withinCandidates)
        fprintf('[clean pipeline] no within-session validation candidates found; skipped\n');
        valWithin = [];
    else
        withinSessionIds = withinCandidates.session_id(1:min(withinValidationMaxSessions, height(withinCandidates)));
        valWithin = ibl_validate_within_session_subset_stitching(qc, cfg, ...
            'SessionIds', withinSessionIds, ...
            'MaxSessions', withinValidationMaxSessions, ...
            'MinAreas', withinValidationMinAreas, ...
            'MaxAreas', withinValidationMaxAreas, ...
            'MinDuration', withinValidationMinDuration, ...
            'NumObservations', withinValidationNumObservations, ...
            'OverlapFraction', withinValidationOverlapFraction, ...
            'NumRepeats', withinValidationNumRepeats, ...
            'Seed', withinValidationSeed, ...
            'Parallel', withinValidationParallel, ...
            'ComputePDC', true, ...
            'MakeFigure', false, ...
            'OutputFile', fullfile(withinDir, 'within_session_subset_stitching_20repeats.mat'));
        writetable(valWithin.summary_table, fullfile(withinDir, 'within_session_subset_stitching_20repeats_summary.csv'));
        local_write_within_validation_summary(valWithin, fullfile(withinDir, 'within_session_subset_stitching_20repeats_session_summary.csv'));
    end
else
    valWithin = [];
    fprintf('[clean pipeline] within-session validation disabled\n');
end

%%
% -------------------------------------------------------------------------
% Step 3: cross-area CSD using the matching Step-1 PC loadings.
% Saves one .mat per session in cfg.cross_spectra_dir. Used by Steps 4-6.
% -------------------------------------------------------------------------
fprintf('\n[clean pipeline] Step 3/7: cross-area CSD with matching PCs\n');
ibl_compute_cross_area_csd(qc, cfg);
local_write_cross_spectra_summary(cfg, fullfile(runRoot, '03_cross_spectra'));

%%
% -------------------------------------------------------------------------
% Step 4: whole-session reproducibility filter.
% Saves filt. Used by Step 5 if enabled.
% -------------------------------------------------------------------------
fprintf('\n[clean pipeline] Step 4/7: session reproducibility filter\n');
filterFile = fullfile(runRoot, '04_session_filter', 'session_repro_filter.mat');
if applySessionReproFilter
    filt = ibl_filter_sessions_by_reproducibility(cfg, ...
        'Band', analysisBand, ...
        'Metric', 'coherence_abs', ...
        'MinSessionsPerPair', reproMinSessionsPerPair, ...
        'MinPairsPerSession', reproMinPairsPerSession, ...
        'NumNull', reproNumNull, ...
        'ThresholdMode', char(reproThresholdMode), ...
        'RobustMadK', reproRobustMadK, ...
        'MakeFigure', true, ...
        'OutputFile', filterFile);
    writetable(filt.session_table, fullfile(runRoot, '04_session_filter', 'session_repro_filter_sessions.csv'));
    cfg.session_filter_file = filterFile;
else
    filt = [];
    cfg.session_filter_file = "";
    fprintf('[clean pipeline] session reproducibility filter disabled\n');
end

% -------------------------------------------------------------------------
% Step 5: derive CSD-correlation threshold, then reliable component.
% Saves plan. Used by Step 6 as the session include list and area universe.
% -------------------------------------------------------------------------
fprintf('\n[clean pipeline] Step 5/7: derive CSD threshold and plan reliable component\n');
planFile = fullfile(runRoot, '05_reliable_plan', 'reliable_plan.mat');
thresholdFile = fullfile(runRoot, '05_reliable_plan', 'csd_threshold_from_validation.mat');
if deriveCsdCorrThreshold
    crossCand = ibl_find_cross_session_validation_candidates(qc, cfg, ...
        'MinTargetAreas', crossValidationMinTargetAreas, ...
        'MaxTargetAreas', crossValidationMaxTargetAreas, ...
        'MinSourceAreas', crossValidationMinSourceAreas, ...
        'MinDuration', crossValidationMinDuration, ...
        'MaxSources', 4, ...
        'MinCoverage', crossValidationMinCoverage, ...
        'MinAnchorAreas', crossValidationMinAnchorAreas, ...
        'MinAnchorFraction', crossValidationMinAnchorFraction, ...
        'MinObservedPairFraction', crossValidationMinObservedPairFraction, ...
        'MaxCompletedPairFraction', crossValidationMaxCompletedPairFraction, ...
        'RequireConnectedSources', crossValidationRequireConnectedSources, ...
        'MinSourceOverlapAreas', crossValidationMinSourceOverlapAreas, ...
        'SessionFilterFile', cfg.session_filter_file, ...
        'TopN', crossValidationTopN, ...
        'Verbose', true);
    if istable(crossCand.candidates)
        writetable(crossCand.candidates, fullfile(runRoot, '05_reliable_plan', 'threshold_validation_candidates.csv'));
    end
    nCandidates = numel(crossCand.sets);
    if nCandidates >= 3
        candidateIdx = 1:min(crossValidationMaxCandidates, nCandidates);
        valThreshold = ibl_validate_cross_session_stitching(crossCand, cfg, ...
            'CandidateIndex', candidateIdx, ...
            'ComputePDC', true, ...
            'MakeFigure', false, ...
            'OutputFile', fullfile(runRoot, '05_reliable_plan', 'threshold_cross_session_validation.mat'));
        thresholdFit = local_fit_csd_threshold(valThreshold.summary_table, thresholdOutcome, csdCorrThresholdFallback);
        csdCorrThreshold = thresholdFit.threshold;
        save(thresholdFile, 'thresholdFit', 'valThreshold', 'crossCand', '-v7.3');
        writetable(thresholdFit.fit_table, fullfile(runRoot, '05_reliable_plan', 'csd_threshold_fit.csv'));
        fprintf('[clean pipeline] derived csdCorrThreshold=%.3f from %s lower 95%% CI crossing 0\n', ...
            csdCorrThreshold, char(thresholdOutcome));
    else
        csdCorrThreshold = csdCorrThresholdFallback;
        thresholdFit = struct('threshold', csdCorrThreshold, ...
            'status', "fallback_too_few_cross_validation_candidates", ...
            'fallback_threshold', csdCorrThresholdFallback);
        save(thresholdFile, 'thresholdFit', 'crossCand', '-v7.3');
        fprintf('[clean pipeline] only %d validation candidates; using fallback csdCorrThreshold=%.3f\n', ...
            nCandidates, csdCorrThreshold);
    end
else
    csdCorrThreshold = csdCorrThresholdFallback;
    thresholdFit = struct('threshold', csdCorrThreshold, ...
        'status', "fixed_fallback_threshold", ...
        'fallback_threshold', csdCorrThresholdFallback);
    save(thresholdFile, 'thresholdFit', '-v7.3');
    fprintf('[clean pipeline] using fixed csdCorrThreshold=%.3f\n', csdCorrThreshold);
end

plan = ibl_plan_reliable_session_growth(cfg, csdCorrThreshold, ...
    'Band', analysisBand, ...
    'Metric', 'coherence_abs', ...
    'MinOverlapAreas', planMinOverlapAreas, ...
    'MinOverlapPairs', planMinOverlapPairs, ...
    'OverlapComponentSelection', char(planOverlapComponentSelection), ...
    'ReliableComponentSelection', char(planReliableComponentSelection), ...
    'SessionFilterFile', cfg.session_filter_file, ...
    'OutputFile', planFile, ...
    'MakeFigure', true);
local_write_plan_summary(plan, fullfile(runRoot, '05_reliable_plan'));
coverageSweep = ibl_sweep_pair_observation_coverage(plan, [0.3 0.4 0.5 0.6], ...
    'OutputFile', fullfile(runRoot, '05_reliable_plan', 'observed_pair_coverage_sweep.mat'), ...
    'Verbose', true);

%%
% -------------------------------------------------------------------------
% Step 6: stitch selected session CSDs and estimate PDC.
% Saves cfg.stitch_file. Used by Step 7.
% -------------------------------------------------------------------------
fprintf('\n[clean pipeline] Step 6/7: stitch selected CSDs and compute PDC\n');
cfg.session_include_ids = plan.selected_session_ids(:);
cfg.stitch_file = fullfile(runRoot, '06_stitch', 'stitched_pdc.mat');
finalQc = qc;
finalQc.largest_component = plan.area_names(:);
finalQc.qualifying_sessions = local_filter_scan_sessions(qc.qualifying_sessions, plan.selected_session_ids);
result = ibl_stitch_saved_spectra(cfg, finalQc);
local_write_stitch_summary(result, fullfile(runRoot, '06_stitch'));

%%
tmp = load(fullfile(runRoot, '06_stitch', 'stitched_pdc.mat'), 'result');
result = tmp.result;

tmpCfg = load(fullfile(runRoot, 'cfg_clean_pipeline.mat'), 'cfg');
cfg = tmpCfg.cfg;

allenMatrixCsv = fullfile(pwd, 'ibl_output', 'allen_tracing', ...
    'allen_tracing_matrix__projection_energy__max.csv');

blockCtrl = ibl_block_pdc_control(result, cfg, allenMatrixCsv, ...
    'Band', [1 80], ...
    'RestrictToAllenSupported', true, ...
    'AllenPositiveThreshold', 500, ...
    'Regularizer', 'eigfloor', ...
    'Lambda', 0, ...
    'MakeFigure', true, ...
    'OutputFile', fullfile(runRoot, '07_infoflow_allen', 'block_pdc_control.mat'));

blockCtrl.summary_table
writetable(blockCtrl.summary_table, ...
    fullfile(runRoot, '07_infoflow_allen', 'block_pdc_control_vs_stitched_observed.csv'));

%%
% -------------------------------------------------------------------------
% Step 7: information flow and Allen comparison.
% Saves infoflow and optional Allen comparison. This is the final figure input.
% -------------------------------------------------------------------------
fprintf('\n[clean pipeline] Step 7/7: info-flow and Allen comparison\n');
fMask = result.freqs >= analysisBand(1) & result.freqs <= min(analysisBand(2), result.freqs(end));
infoflow = ibl_pdc_to_infoflow(result.PDC(:,:,fMask));
observationMask = local_observation_mask(result);
save(fullfile(runRoot, '07_infoflow_allen', 'infoflow.mat'), ...
    'infoflow', 'observationMask', 'analysisBand', 'result', '-v7.3');
local_write_matrix_csv(infoflow, string(result.area_names(:)), ...
    fullfile(runRoot, '07_infoflow_allen', 'infoflow.csv'));
local_write_matrix_csv(double(observationMask), string(result.area_names(:)), ...
    fullfile(runRoot, '07_infoflow_allen', 'direct_coobservation_mask.csv'));

allenMatrixCsv = fullfile(pwd, 'ibl_output', 'allen_tracing', ...
    'allen_tracing_matrix__projection_density__max.csv');

if exist(allenMatrixCsv, 'file') == 2
    cmp = ibl_compare_allen_tracing(result, allenMatrixCsv, ...
        'Band', analysisBand, ...
        'ObservationMask', observationMask, ...
        'MakeMatrixFigure', true, ...
        'MakeRocFigure', true);
    save(fullfile(runRoot, '07_infoflow_allen', 'allen_comparison.mat'), 'cmp', '-v7.3');
    local_write_allen_summary(cmp, fullfile(runRoot, '07_infoflow_allen', 'allen_comparison_summary.csv'));
    blockCtrl = ibl_block_pdc_control(result, cfg, allenMatrixCsv, ...
        'Band', analysisBand, ...
        'RestrictToAllenSupported', true, ...
        'AllenPositiveThreshold', 0, ...
        'Regularizer', 'eigfloor', ...
        'Lambda', 0, ...
        'MakeFigure', true, ...
        'OutputFile', fullfile(runRoot, '07_infoflow_allen', 'block_pdc_control.mat'));
    writetable(blockCtrl.summary_table, ...
        fullfile(runRoot, '07_infoflow_allen', 'block_pdc_control_vs_stitched_observed.csv'));
    writetable(blockCtrl.block_table, ...
        fullfile(runRoot, '07_infoflow_allen', 'block_pdc_control_sessions.csv'));
    local_write_matrix_csv(blockCtrl.infoflow_observed_mean, string(result.area_names(:)), ...
        fullfile(runRoot, '07_infoflow_allen', 'block_pdc_control_infoflow_observed_mean.csv'));
else
    cmp = [];
    blockCtrl = [];
    fprintf('[clean pipeline] Allen matrix not found, skipped: %s\n', char(allenMatrixCsv));
end

save(fullfile(runRoot, 'clean_pipeline_workspace.mat'), ...
    'cfg', 'scan', 'qc', 'valWithin', 'filt', 'plan', 'coverageSweep', 'result', 'infoflow', ...
    'observationMask', 'cmp', 'blockCtrl', '-v7.3');

fprintf('\n[clean pipeline] complete: %s\n', runRoot);

function cfg = local_apply_clean_paths(cfg, runRoot)
cfg.scan_file = fullfile(runRoot, '00_scan', 'session_scan.mat');
cfg.area_pca_dir = fullfile(runRoot, '01_area_pca', 'area_pca');
cfg.pca_qc_file = fullfile(runRoot, '02_pca_qc', 'pca_qc.mat');
cfg.cross_spectra_dir = fullfile(runRoot, '03_cross_spectra', 'cross_spectra');
cfg.session_spectra_dir = cfg.cross_spectra_dir;
cfg.stitch_file = fullfile(runRoot, '06_stitch', 'stitched_pdc.mat');
cfg.session_filter_file = "";
end

function local_make_clean_dirs(cfg)
dirs = { ...
    fileparts(cfg.scan_file), cfg.area_pca_dir, fileparts(cfg.pca_qc_file), ...
    cfg.cross_spectra_dir, ...
    fullfile(cfg.output_root, '04_session_filter'), ...
    fullfile(cfg.output_root, '02_within_session_validation'), ...
    fullfile(cfg.output_root, '05_reliable_plan'), ...
    fullfile(cfg.output_root, '06_stitch'), ...
    fullfile(cfg.output_root, '07_infoflow_allen')};
for i = 1:numel(dirs)
    if exist(dirs{i}, 'dir') ~= 7
        mkdir(dirs{i});
    end
end
end

function local_write_config_summary(cfg, runRoot, runName, analysisBand, applyFilter, ...
    minSessPair, minPairsSess, threshMode, robustK, numNull, deriveCsdThresh, ...
    csdThreshFallback, thresholdOutcome, ...
    minOverlapAreas, minOverlapPairs, overlapSelection, reliableSelection, allenCsv)
fid = fopen(fullfile(runRoot, 'pipeline_config.txt'), 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'run_name: %s\n', char(runName));
fprintf(fid, 'created_at: %s\n', datestr(now));
fprintf(fid, 'cache_root: %s\n', cfg.cache_root);
fprintf(fid, 'output_root: %s\n\n', cfg.output_root);
fprintf(fid, '[neuron/session gates]\n');
fprintf(fid, 'min_sp_dur: %.6g\n', cfg.min_sp_dur);
fprintf(fid, 'min_mean_fr: %.6g\n', cfg.min_mean_fr);
fprintf(fid, 'min_neurons_per_area: %d\n', cfg.min_neurons_per_area);
fprintf(fid, 'max_neurons_per_area: %d\n', cfg.max_neurons_per_area);
fprintf(fid, 'min_pc1_explained: %.6g\n', cfg.min_pc1_explained);
fprintf(fid, 'min_sessions_per_area: %d\n\n', cfg.min_sessions_per_area);
fprintf(fid, '[spectral PCA]\n');
fprintf(fid, 'pc_band: [%g %g]\n', cfg.pc_band(1), cfg.pc_band(2));
fprintf(fid, 'pc_sign_convention: largest_abs_loading_real_positive\n');
fprintf(fid, 'target_n_freqs: %d\n', cfg.target_n_freqs);
fprintf(fid, 'chronux.win: %.6g\n', cfg.chronux.win);
fprintf(fid, 'chronux.tapers: [%g %g]\n\n', cfg.chronux.tapers(1), cfg.chronux.tapers(2));
fprintf(fid, '[within-session validation]\n');
if isfield(cfg, 'within_session_validation')
    wv = cfg.within_session_validation;
    fprintf(fid, 'enabled: %d\n', wv.enabled);
    fprintf(fid, 'num_repeats: %d\n', wv.num_repeats);
    fprintf(fid, 'max_sessions: %d\n', wv.max_sessions);
    fprintf(fid, 'min_areas: %d\n', wv.min_areas);
    fprintf(fid, 'max_areas: %d\n', wv.max_areas);
    fprintf(fid, 'min_duration: %.6g\n', wv.min_duration);
    fprintf(fid, 'num_observations: %d\n', wv.num_observations);
    fprintf(fid, 'overlap_fraction: %.6g\n', wv.overlap_fraction);
    fprintf(fid, 'seed: %d\n', wv.seed);
    fprintf(fid, 'parallel: %d\n\n', wv.parallel);
end
fprintf(fid, '[session reproducibility]\n');
fprintf(fid, 'enabled: %d\n', applyFilter);
fprintf(fid, 'analysis_band: [%g %g]\n', analysisBand(1), analysisBand(2));
fprintf(fid, 'min_sessions_per_pair: %d\n', minSessPair);
fprintf(fid, 'min_pairs_per_session: %d\n', minPairsSess);
fprintf(fid, 'threshold_mode: %s\n', char(threshMode));
fprintf(fid, 'robust_mad_k: %.6g\n', robustK);
fprintf(fid, 'num_null: %d\n\n', numNull);
fprintf(fid, '[reliable component]\n');
fprintf(fid, 'derive_csd_corr_threshold: %d\n', deriveCsdThresh);
fprintf(fid, 'csd_corr_threshold_fallback: %.6g\n', csdThreshFallback);
fprintf(fid, 'threshold_outcome: %s\n', char(thresholdOutcome));
if isfield(cfg, 'threshold_validation')
    tv = cfg.threshold_validation;
    fprintf(fid, 'threshold_validation_top_n: %d\n', tv.top_n);
    fprintf(fid, 'threshold_validation_max_candidates: %d\n', tv.max_candidates);
    fprintf(fid, 'threshold_validation_min_target_areas: %d\n', tv.min_target_areas);
    fprintf(fid, 'threshold_validation_max_target_areas: %d\n', tv.max_target_areas);
    fprintf(fid, 'threshold_validation_min_source_areas: %d\n', tv.min_source_areas);
    fprintf(fid, 'threshold_validation_min_duration: %.6g\n', tv.min_duration);
    fprintf(fid, 'threshold_validation_min_coverage: %.6g\n', tv.min_coverage);
    fprintf(fid, 'threshold_validation_min_anchor_areas: %d\n', tv.min_anchor_areas);
    fprintf(fid, 'threshold_validation_min_anchor_fraction: %.6g\n', tv.min_anchor_fraction);
    fprintf(fid, 'threshold_validation_min_observed_pair_fraction: %.6g\n', tv.min_observed_pair_fraction);
    fprintf(fid, 'threshold_validation_max_completed_pair_fraction: %.6g\n', tv.max_completed_pair_fraction);
    fprintf(fid, 'threshold_validation_require_connected_sources: %d\n', tv.require_connected_sources);
    fprintf(fid, 'threshold_validation_min_source_overlap_areas: %d\n', tv.min_source_overlap_areas);
end
fprintf(fid, 'min_overlap_areas: %d\n', minOverlapAreas);
fprintf(fid, 'min_overlap_pairs: %d\n', minOverlapPairs);
fprintf(fid, 'overlap_component_selection: %s\n', char(overlapSelection));
fprintf(fid, 'reliable_component_selection: %s\n\n', char(reliableSelection));
fprintf(fid, '[stitching]\n');
fprintf(fid, 'normalize: %s\n', cfg.stitch.normalize);
fprintf(fid, 'regularizer: %s\n', cfg.stitch.regularizer);
fprintf(fid, 'lambda: %.6g\n', cfg.stitch.lambda);
fprintf(fid, 'parallel: %d\n\n', cfg.stitch.parallel);
fprintf(fid, '[allen]\n');
fprintf(fid, 'allen_matrix_csv: %s\n', char(allenCsv));
save(fullfile(runRoot, 'cfg_clean_pipeline.mat'), 'cfg');
end

function local_write_scan_summary(scan, outDir)
sessions = scan.sessions(:);
T = local_session_struct_table(sessions);
writetable(T, fullfile(outDir, 'scan_all_sessions.csv'));
Tq = local_session_struct_table(scan.qualifying_sessions(:));
writetable(Tq, fullfile(outDir, 'scan_qualifying_sessions.csv'));
area = scan.area_names(:);
in_largest_component = ismember(area, scan.largest_component(:));
writetable(table(area, in_largest_component), fullfile(outDir, 'scan_area_component.csv'));
end

function T = local_session_struct_table(sessions)
n = numel(sessions);
session_id = strings(n,1); lab = strings(n,1); subject = strings(n,1);
date = strings(n,1); number = strings(n,1); duration_s = nan(n,1);
n_areas = nan(n,1); areas = strings(n,1);
for i = 1:n
    session_id(i) = string(sessions(i).session_id);
    [lab0, subject0, date0, number0] = local_parse_session_id(session_id(i));
    lab(i) = local_get_session_string(sessions(i), 'lab', lab0);
    subject(i) = local_get_session_string(sessions(i), 'subject', subject0);
    date(i) = local_get_session_string(sessions(i), 'date', date0);
    number(i) = local_get_session_string(sessions(i), 'number', number0);
    duration_s(i) = local_get_session_double(sessions(i), 'sp_dur', NaN);
    n_areas(i) = local_get_session_double(sessions(i), 'n_areas', NaN);
    if isfield(sessions(i), 'areas')
        areas(i) = strjoin(string(sessions(i).areas(:)).', "; ");
    elseif isfield(sessions(i), 'area_names')
        areas(i) = strjoin(string(sessions(i).area_names(:)).', "; ");
    else
        areas(i) = "";
    end
end
T = table(session_id, lab, subject, date, number, duration_s, n_areas, areas);
end

function [lab, subject, date, number] = local_parse_session_id(sessionId)
parts = split(string(sessionId), "__");
lab = ""; subject = ""; date = ""; number = "";
if numel(parts) >= 1, lab = parts(1); end
if numel(parts) >= 2, subject = parts(2); end
if numel(parts) >= 3, date = parts(3); end
if numel(parts) >= 4, number = parts(4); end
end

function value = local_get_session_string(s, fieldName, defaultValue)
if isfield(s, fieldName)
    value = string(s.(fieldName));
else
    value = string(defaultValue);
end
end

function value = local_get_session_double(s, fieldName, defaultValue)
if isfield(s, fieldName)
    value = double(s.(fieldName));
else
    value = defaultValue;
end
end

function local_write_area_pca_summary(cfg, outDir)
files = dir(fullfile(cfg.area_pca_dir, '*.mat'));
rows = struct([]);
for k = 1:numel(files)
    tmp = load(fullfile(files(k).folder, files(k).name), 'summary');
    s = tmp.summary;
    names = string(s.area_names(:));
    for a = 1:numel(names)
        row.session_id = string(s.session_id);
        row.area = names(a);
        row.n_neurons = double(s.n_neurons(a));
        row.duration_s = double(s.sp_dur);
        row.mean_pc1_explained = double(s.mean_pc1_explained(a));
        row.pass_qc_area = logical(s.pass_qc_area(a));
        row.qc_reason = string(s.qc_reason(a));
        row.pc_sign_convention = string(s.pc_sign_convention);
        rows = local_append_struct(rows, row);
    end
end
writetable(local_struct_to_table(rows), fullfile(outDir, 'area_pca_summary.csv'));
end

function local_write_qc_summary(qc, outDir)
area = qc.area_names(:);
in_largest_component = ismember(area, qc.largest_component(:));
writetable(table(area, in_largest_component), fullfile(outDir, 'pca_qc_areas.csv'));
writetable(local_session_struct_table(qc.qualifying_sessions(:)), fullfile(outDir, 'pca_qc_qualifying_sessions.csv'));
end

function local_write_within_validation_summary(val, outFile)
if isempty(val) || ~isfield(val, 'summary_table') || isempty(val.summary_table)
    writetable(table(), outFile);
    return
end
T = val.summary_table;
session_id = unique(T.session_id, 'stable');
n = numel(session_id);
n_repeats = zeros(n,1);
info_observed_corr_median = nan(n,1);
info_completed_corr_median = nan(n,1);
info_completed_corr_iqr_low = nan(n,1);
info_completed_corr_iqr_high = nan(n,1);
info_observed_rel_rmse_median = nan(n,1);
info_completed_rel_rmse_median = nan(n,1);
for i = 1:n
    idx = T.session_id == session_id(i);
    n_repeats(i) = nnz(idx);
    if ismember('info_observed_corr', T.Properties.VariableNames)
        info_observed_corr_median(i) = median(T.info_observed_corr(idx), 'omitnan');
        info_completed_corr_median(i) = median(T.info_completed_corr(idx), 'omitnan');
        info_completed_corr_iqr_low(i) = prctile(T.info_completed_corr(idx), 25);
        info_completed_corr_iqr_high(i) = prctile(T.info_completed_corr(idx), 75);
        info_observed_rel_rmse_median(i) = median(T.info_observed_rel_rmse(idx), 'omitnan');
        info_completed_rel_rmse_median(i) = median(T.info_completed_rel_rmse(idx), 'omitnan');
    end
end
summaryTable = table(session_id, n_repeats, info_observed_corr_median, ...
    info_completed_corr_median, info_completed_corr_iqr_low, ...
    info_completed_corr_iqr_high, info_observed_rel_rmse_median, ...
    info_completed_rel_rmse_median);
writetable(summaryTable, outFile);
end

function local_write_cross_spectra_summary(cfg, outDir)
files = dir(fullfile(cfg.cross_spectra_dir, '*.mat'));
session_id = strings(numel(files),1);
n_areas = zeros(numel(files),1);
n_freqs = zeros(numel(files),1);
areas = strings(numel(files),1);
for k = 1:numel(files)
    tmp = load(fullfile(files(k).folder, files(k).name), 'summary');
    s = tmp.summary;
    session_id(k) = string(s.session_id);
    n_areas(k) = numel(s.area_names);
    n_freqs(k) = numel(s.freqs);
    areas(k) = strjoin(string(s.area_names(:)).', "; ");
end
writetable(table(session_id, n_areas, n_freqs, areas), fullfile(outDir, 'cross_spectra_summary.csv'));
end

function local_write_plan_summary(plan, outDir)
writetable(plan.selected_table, fullfile(outDir, 'selected_sessions.csv'));
writetable(plan.eligible_session_table, fullfile(outDir, 'eligible_sessions.csv'));
writetable(plan.reliable_component_table, fullfile(outDir, 'reliable_components.csv'));
area = plan.area_names(:);
support_count = plan.area_support_count(:);
writetable(table(area, support_count), fullfile(outDir, 'selected_area_support.csv'));
end

function thresholdFit = local_fit_csd_threshold(T, outcomeName, fallbackThreshold)
x = T.csd_observed_abs_corr;
if ~ismember(outcomeName, string(T.Properties.VariableNames))
    error('Threshold outcome %s is not in validation summary table.', char(outcomeName));
end
y = T.(char(outcomeName));
good = isfinite(x) & isfinite(y);
if nnz(good) < 3 || range(x(good)) == 0
    thresholdFit.threshold = fallbackThreshold;
    thresholdFit.status = "fallback_too_few_points";
    thresholdFit.fallback_threshold = fallbackThreshold;
    thresholdFit.n_points = nnz(good);
    thresholdFit.fit_table = table();
    return
end

mdl = fitlm(x(good), y(good));
xGrid = linspace(max(0, min(x(good))), min(1, max(x(good))), 400).';
[yHat, yCI] = predict(mdl, xGrid, 'Alpha', 0.05, 'Prediction', 'curve');
lower95 = yCI(:,1);
idx = find(lower95 >= 0, 1, 'first');
if isempty(idx)
    threshold = fallbackThreshold;
    status = "fallback_lower95_never_crosses_zero";
elseif idx == 1
    threshold = xGrid(1);
    status = "derived_lower95_already_nonnegative";
else
    x0 = xGrid(idx-1);
    x1 = xGrid(idx);
    y0 = lower95(idx-1);
    y1 = lower95(idx);
    threshold = x0 + (0 - y0) * (x1 - x0) / max(y1 - y0, eps);
    status = "derived_lower95_crosses_zero";
end
threshold = max(0, min(1, threshold));

thresholdFit.threshold = threshold;
thresholdFit.status = status;
thresholdFit.fallback_threshold = fallbackThreshold;
thresholdFit.outcome = string(outcomeName);
thresholdFit.n_points = nnz(good);
thresholdFit.model = mdl;
thresholdFit.fit_table = table(xGrid, yHat, lower95, yCI(:,2), ...
    'VariableNames', {'csd_observed_abs_corr', 'predicted_infoflow_corr', ...
    'lower95', 'upper95'});
end

function sessionsOut = local_filter_scan_sessions(sessionsIn, ids)
ids = string(ids(:));
keep = false(numel(sessionsIn), 1);
for i = 1:numel(sessionsIn)
    keep(i) = any(string(sessionsIn(i).session_id) == ids);
end
sessionsOut = sessionsIn(keep);
end

function local_write_stitch_summary(result, outDir)
area = string(result.area_names(:));
support_count = zeros(numel(area),1);
for s = 1:numel(result.kept_sessions)
    support_count = support_count + ismember(area, string(result.kept_sessions(s).areas(:)));
end
writetable(table(area, support_count), fullfile(outDir, 'stitched_area_support.csv'));
session_id = strings(numel(result.kept_sessions),1);
n_areas = zeros(numel(result.kept_sessions),1);
areas = strings(numel(result.kept_sessions),1);
for s = 1:numel(result.kept_sessions)
    session_id(s) = string(result.kept_sessions(s).session_id);
    n_areas(s) = numel(result.kept_sessions(s).areas);
    areas(s) = strjoin(string(result.kept_sessions(s).areas(:)).', "; ");
end
writetable(table(session_id, n_areas, areas), fullfile(outDir, 'stitched_sessions.csv'));
end

function obs = local_observation_mask(result)
n = numel(result.area_names);
obs = false(n);
for r = 1:numel(result.recset)
    idx = result.recset{r};
    obs(idx, idx) = true;
end
obs(1:n+1:end) = false;
end

function local_write_matrix_csv(M, labels, outFile)
T = array2table(M, 'VariableNames', matlab.lang.makeValidName(cellstr(labels)));
T = addvars(T, labels(:), 'Before', 1, 'NewVariableNames', 'target_area');
writetable(T, outFile);
end

function local_write_allen_summary(cmp, outFile)
metric = ["roc_auc"; "observed_roc_auc"; "never_observed_roc_auc"];
value = [local_get_field(cmp, 'roc_auc'); local_get_field(cmp, 'observed_roc_auc'); ...
    local_get_field(cmp, 'never_observed_roc_auc')];
writetable(table(metric, value), outFile);
end

function x = local_get_field(S, name)
if isfield(S, name)
    x = S.(name);
else
    x = nan;
end
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
