function cfg = ibl_default_config(cacheRoot, outputRoot)
%IBL_DEFAULT_CONFIG Default configuration for the IBL spontaneous pipeline.
%
%   cfg = ibl_default_config(cacheRoot, outputRoot)

if nargin < 1 || isempty(cacheRoot)
    cacheRoot = '/Volumes/Extreme SSD/data/neuropixel/ONE/openalyx.internationalbrainlab.org';
end
if nargin < 2 || isempty(outputRoot)
    outputRoot = fullfile(pwd, 'ibl_output');
end

cfg.cache_root = cacheRoot;
cfg.output_root = outputRoot;
cfg.session_spectra_dir = fullfile(outputRoot, 'session_spectra'); % legacy
cfg.area_pca_dir = fullfile(outputRoot, 'area_pca');
cfg.cross_spectra_dir = fullfile(outputRoot, 'cross_spectra');
cfg.scan_file = fullfile(outputRoot, 'session_scan.mat');
cfg.pca_qc_file = fullfile(outputRoot, 'pca_qc.mat');
cfg.stitch_file = fullfile(outputRoot, 'stitched_pdc.mat');
cfg.session_filter_file = "";  % optional output from ibl_filter_sessions_by_reproducibility
cfg.verbose = true;
cfg.progress_every = 10;
cfg.max_sessions = inf;

cfg.min_sp_dur = 300;           % seconds
cfg.min_neurons_per_area = 5;
cfg.min_mean_fr = 0.2;          % Hz
cfg.max_neurons_per_area = 30;
cfg.n_pcs = 3;
cfg.min_pc1_explained = 0.10;   % area/session quality threshold
cfg.pc_band = [1 100];          % Hz band used to define fixed area loading
cfg.target_n_freqs = 128;       % reduce frequency grid before downstream analysis
cfg.min_sessions_per_area = 3;  % post-PCA-QC support threshold

cfg.chronux.Fs = 1000;
cfg.chronux.tapers = [4 7];
cfg.chronux.pad = 0;
cfg.chronux.fpass = [0 100];
cfg.chronux.err = [0 0];
cfg.chronux.trialave = 0;
cfg.chronux.win = 4.0;          % seconds
cfg.chronux.maxiter = 100;
cfg.chronux.tol = 1e-3;

cfg.stitch.method = 'nnm';
cfg.stitch.normalize = 'none';  % 'none' or 'coherence'
cfg.stitch.lambda = 0;
cfg.stitch.regularizer = 'eigfloor';  % fast SPD repair before Wilson/PDC
cfg.stitch.min_eig_rel = 1e-6;
cfg.stitch.min_eig_abs = 1e-8;
cfg.stitch.glasso_fallback_lambda = 0.01;
cfg.stitch.glasso_maxiter = 30;
cfg.stitch.glasso_tol = 1e-4;
cfg.stitch.glasso_retry_factor = 5;
cfg.stitch.glasso_standardize = true;
cfg.stitch.parallel = false;
cfg.stitch.glasso_diagnostics = false;
cfg.stitch.glasso_live_print = false;

cfg.area_map = ibl_area_map();

if ~exist(cfg.output_root, 'dir')
    mkdir(cfg.output_root);
end
if ~exist(cfg.session_spectra_dir, 'dir')
    mkdir(cfg.session_spectra_dir);
end
if ~exist(cfg.area_pca_dir, 'dir')
    mkdir(cfg.area_pca_dir);
end
if ~exist(cfg.cross_spectra_dir, 'dir')
    mkdir(cfg.cross_spectra_dir);
end
end
