function summaries = ibl_compute_area_pca(scan, cfg)
%IBL_COMPUTE_AREA_PCA Compute within-area spectra and spectral PCA only.

assert(exist('CrossSpecMatpt', 'file') == 2, ...
    'Chronux CrossSpecMatpt not found. Add Chronux to the MATLAB path.');

sessions = scan.qualifying_sessions;
regionTable = scan.region_table;
allowedAreas = scan.largest_component;
if isfinite(cfg.max_sessions)
    sessions = sessions(1:min(numel(sessions), cfg.max_sessions));
end

summaries = struct('session_id', {}, 'file', {}, 'area_names', {});
if cfg.verbose
    fprintf('[ibl_compute_area_pca] Processing %d qualifying sessions\n', numel(sessions));
end
runTimer = tic;

for k = 1:numel(sessions)
    sessionId = sessions(k).session_id;
    outFile = fullfile(cfg.area_pca_dir, [sessionId '.mat']);
    if cfg.verbose
        fprintf('[ibl_compute_area_pca] %4d/%4d  %s\n', k, numel(sessions), sessionId);
    end
    if exist(outFile, 'file') == 2
        tmp = load(outFile, 'summary');
        rec.session_id = tmp.summary.session_id;
        rec.file = outFile;
        rec.area_names = string(tmp.summary.area_names);
        summaries(end+1) = rec; %#ok<AGROW>
        if cfg.verbose
            fprintf('  cached: %s\n', outFile);
        end
        continue;
    end

    areaStruct = ibl_load_area_spikes(sessions(k).session_root, regionTable, cfg);
    keep = ismember(string({areaStruct.area}), allowedAreas);
    areaStruct = areaStruct(keep);
    if numel(areaStruct) < 2
        if cfg.verbose
            fprintf('  skip: only %d allowed areas after reload\n', numel(areaStruct));
        end
        continue;
    end

    areaNames = string({areaStruct.area});
    [areaNames, ord] = sort(areaNames);
    areaStruct = areaStruct(ord);
    nAreas = numel(areaStruct);
    Tmax = sessions(k).sp_dur;
    if cfg.verbose
        fprintf('  areas=%d  duration=%.1fs\n', nAreas, Tmax);
    end

    eigvals = cell(nAreas,1);
    eigvecs = cell(nAreas,1);
    explVar = cell(nAreas,1);
    fixedLoadings = cell(nAreas,1);
    fixedEigvals = zeros(nAreas, 1);
    autoSpectrum = [];
    freqs = [];
    passQCArea = false(nAreas,1);
    qcReason = strings(nAreas,1);

    for a = 1:nAreas
        data = local_times_to_struct(areaStruct(a).spike_times);
        [Saa, f] = local_point_spectrum(data, Tmax, cfg.chronux);
        [Saa, f] = ibl_reduce_frequency_grid(Saa, f, cfg.target_n_freqs);
        if isempty(freqs)
            freqs = f;
            nFreq = numel(f);
            autoSpectrum = zeros(nAreas, nFreq);
        end
        [ev, evec, expl] = ibl_spectral_pca(Saa, cfg.n_pcs);
        [vfix, dfix] = ibl_fixed_spectral_loading(Saa, freqs, cfg.pc_band);
        eigvals{a} = ev;
        eigvecs{a} = evec;
        explVar{a} = expl;
        fixedLoadings{a} = vfix;
        fixedEigvals(a) = dfix;
        autoSpectrum(a,:) = ev(1,:);

        pc1score = mean(expl(1, freqs >= cfg.pc_band(1) & freqs <= min(cfg.pc_band(2), freqs(end))), 'omitnan');
        if ~isfinite(pc1score)
            passQCArea(a) = false;
            qcReason(a) = "nonfinite_pc1";
        elseif pc1score < cfg.min_pc1_explained
            passQCArea(a) = false;
            qcReason(a) = "low_pc1";
        else
            passQCArea(a) = true;
            qcReason(a) = "pass";
        end

        if cfg.verbose
            fprintf('    area %2d/%2d  %-12s  neurons=%2d  pc1=%.2f  qc=%s\n', ...
                a, nAreas, char(areaNames(a)), areaStruct(a).n_neurons, pc1score, char(qcReason(a)));
        end
    end

    summary.session_id = sessionId;
    summary.session_root = sessions(k).session_root;
    summary.area_names = areaNames;
    summary.n_neurons = [areaStruct.n_neurons]';
    summary.firing_rates = {areaStruct.firing_rates};
    summary.freqs = freqs;
    summary.auto_spectrum = autoSpectrum;
    summary.eigenvalues = eigvals;
    summary.eigenvectors = eigvecs;
    summary.expl_var = explVar;
    summary.fixed_loadings = fixedLoadings;
    summary.fixed_loading_eigvals = fixedEigvals;
    summary.pc_sign_convention = "largest_abs_loading_real_positive";
    summary.pc_band = cfg.pc_band;
    summary.mean_pc1_explained = local_mean_pc1(explVar, freqs, cfg.pc_band);
    summary.pass_qc_area = passQCArea;
    summary.qc_reason = qcReason;
    summary.sp_dur = Tmax;
    summary.chronux = cfg.chronux;
    summary.target_n_freqs = cfg.target_n_freqs;

    save(outFile, 'summary', '-v7.3');
    if cfg.verbose
        fprintf('  saved: %s\n', outFile);
        fprintf('  elapsed total: %.1fs\n', toc(runTimer));
    end
    rec.session_id = sessionId;
    rec.file = outFile;
    rec.area_names = areaNames;
    summaries(end+1) = rec; %#ok<AGROW>
end
end


function data = local_times_to_struct(timesCell)
n = numel(timesCell);
data = struct('times', cell(n,1));
for i = 1:n
    data(i).times = timesCell{i}(:);
end
end


function [S, f] = local_point_spectrum(data, Tmax, chronux)
[pxy, ~, ~, ~, ~, f] = CrossSpecMatpt(data, chronux.win, Tmax, chronux);
S = permute(pxy, [3, 2, 1]);
for i = 1:size(S,3)
    S(:,:,i) = (S(:,:,i) + S(:,:,i)') / 2;
end
S(~isfinite(S)) = 0;
end


function scores = local_mean_pc1(explVar, freqs, band)
n = numel(explVar);
scores = nan(n,1);
mask = freqs >= band(1) & freqs <= min(band(2), freqs(end));
for i = 1:n
    scores(i) = mean(explVar{i}(1, mask), 'omitnan');
end
end
