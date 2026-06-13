function summaries = ibl_compute_cross_area_csd(qc, cfg)
%IBL_COMPUTE_CROSS_AREA_CSD Compute cross-area CSD only for PCA-QC-passed component.

assert(exist('CrossSpecMatpt', 'file') == 2, ...
    'Chronux CrossSpecMatpt not found. Add Chronux to the MATLAB path.');

sessions = qc.qualifying_sessions;
regionTable = qc.scan.region_table;
allowedAreas = qc.largest_component;

summaries = struct('session_id', {}, 'file', {}, 'area_names', {});
if cfg.verbose
    fprintf('[ibl_compute_cross_area_csd] Processing %d PCA-QC sessions\n', numel(sessions));
end

for k = 1:numel(sessions)
    sessionId = sessions(k).session_id;
    outFile = fullfile(cfg.cross_spectra_dir, [sessionId '.mat']);
    if cfg.verbose
        fprintf('[ibl_compute_cross_area_csd] %4d/%4d  %s\n', k, numel(sessions), sessionId);
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

    pcaFile = fullfile(cfg.area_pca_dir, [sessionId '.mat']);
    pcaTmp = load(pcaFile, 'summary');
    pcaSummary = pcaTmp.summary;

    areaStruct = ibl_load_area_spikes(pcaSummary.session_root, regionTable, cfg);
    [allNames, ord] = sort(string({areaStruct.area}));
    areaStruct = areaStruct(ord);

    keep = pcaSummary.pass_qc_area(:) & ismember(string(pcaSummary.area_names(:)), allowedAreas);
    areaNames = string(pcaSummary.area_names(keep));
    if numel(areaNames) < 2
        continue;
    end

    [tf, idxStruct] = ismember(areaNames, allNames);
    areaNames = areaNames(tf);
    idxStruct = idxStruct(tf);
    fixedLoadings = pcaSummary.fixed_loadings(keep);
    fixedLoadings = fixedLoadings(tf);
    if numel(areaNames) < 2
        continue;
    end

    Tmax = pcaSummary.sp_dur;
    freqs = pcaSummary.freqs(:);
    nAreas = numel(areaNames);
    crossSpectrum = zeros(nAreas, nAreas, numel(freqs));
    autoSpectrum = pcaSummary.auto_spectrum(keep, :);
    autoSpectrum = autoSpectrum(tf, :);
    for a = 1:nAreas
        crossSpectrum(a,a,:) = autoSpectrum(a,:);
    end

    for a = 1:nAreas
        for b = (a+1):nAreas
            pairData = [local_times_to_struct(areaStruct(idxStruct(a)).spike_times); ...
                        local_times_to_struct(areaStruct(idxStruct(b)).spike_times)];
            [SabFull, f] = local_point_spectrum(pairData, Tmax, cfg.chronux);
            [SabFull, f] = ibl_reduce_frequency_grid(SabFull, f, cfg.target_n_freqs);
            assert(numel(f) == numel(freqs) && max(abs(f(:)-freqs(:))) < 1e-10, ...
                'Frequency mismatch in cross-area CSD for %s', sessionId);
            na = areaStruct(idxStruct(a)).n_neurons;
            nb = areaStruct(idxStruct(b)).n_neurons;
            Sab = SabFull(1:na, na+1:na+nb, :);
            vab = ibl_project_cross_spectrum_fixed(Sab, fixedLoadings{a}, fixedLoadings{b});
            crossSpectrum(a,b,:) = vab;
            crossSpectrum(b,a,:) = conj(vab);
            if cfg.verbose
                fprintf('    cross (%s, %s)\n', char(areaNames(a)), char(areaNames(b)));
            end
        end
    end

    summary.session_id = sessionId;
    summary.area_names = areaNames;
    summary.freqs = freqs;
    summary.cross_spectrum = crossSpectrum;
    summary.auto_spectrum = autoSpectrum;
    summary.pc_band = pcaSummary.pc_band;
    summary.mean_pc1_explained = pcaSummary.mean_pc1_explained(keep);
    summary.mean_pc1_explained = summary.mean_pc1_explained(tf);
    save(outFile, 'summary', '-v7.3');

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
