function [P, S, Sreg, f, meta] = stitch_spectra_blocks(Sblocks, recset, f, params)
%STITCH_SPECTRA_BLOCKS Stitch precomputed block spectra into a global spectrum.

arguments
    Sblocks (1,:) cell
    recset (1,:) cell
    f (:,1) double
    params struct
end

if ~isfield(params, 'method'); params.method = 'none'; end
if ~isfield(params, 'normalize'); params.normalize = 'none'; end
if ~isfield(params, 'lambda'); params.lambda = 0; end
if ~isfield(params, 'regularizer'); params.regularizer = 'eigfloor'; end
if ~isfield(params, 'min_eig_rel'); params.min_eig_rel = 1e-6; end
if ~isfield(params, 'min_eig_abs'); params.min_eig_abs = 1e-8; end
if ~isfield(params, 'glasso_fallback_lambda'); params.glasso_fallback_lambda = 0.01; end
if ~isfield(params, 'glasso_maxiter'); params.glasso_maxiter = 30; end
if ~isfield(params, 'glasso_tol'); params.glasso_tol = 1e-4; end
if ~isfield(params, 'glasso_retry_factor'); params.glasso_retry_factor = 5; end
if ~isfield(params, 'glasso_standardize'); params.glasso_standardize = true; end
if ~isfield(params, 'parallel'); params.parallel = false; end
if ~isfield(params, 'glasso_diagnostics'); params.glasso_diagnostics = false; end
if ~isfield(params, 'glasso_live_print'); params.glasso_live_print = false; end
if ~isfield(params, 'verbose'); params.verbose = false; end

countMat = meacount_mat(recset);
rawCountMat = countMat;
missMat = (rawCountMat == 0);
K = size(rawCountMat, 1);
nf = numel(f);

if params.verbose
    if isscalar(params.lambda)
        lambdaText = sprintf('%.4g', params.lambda);
    else
        lambdaText = sprintf('[%.4g %.4g]', min(params.lambda(:)), max(params.lambda(:)));
    end
    fprintf('[stitch_spectra_blocks] blocks=%d  nodes=%d  freqs=%d  method=%s  normalize=%s  lambda=%s  reg=%s\n', ...
        numel(Sblocks), K, nf, params.method, params.normalize, lambdaText, params.regularizer);
end

S = zeros(K, K, nf);
autoPower = zeros(K, nf);
autoCount = zeros(K, 1);
for u = 1:numel(Sblocks)
    idx = recset{u};
    Su = Sblocks{u};
    assert(size(Su,1) == numel(idx) && size(Su,2) == numel(idx), ...
        'Block %d size does not match recset.', u);
    assert(size(Su,3) == nf, 'Frequency dimension mismatch for block %d.', u);
    SuRaw = Su;
    Su = local_normalize_block(Su, params.normalize);
    S(idx, idx, :) = S(idx, idx, :) + Su;
    for a = 1:numel(idx)
        autoPower(idx(a), :) = autoPower(idx(a), :) + squeeze(real(SuRaw(a,a,:))).';
        autoCount(idx(a)) = autoCount(idx(a)) + 1;
    end
    if params.verbose
        fprintf('[stitch_spectra_blocks] added block %3d/%3d  size=%d\n', ...
            u, numel(Sblocks), numel(idx));
    end
end

countMatSafe = rawCountMat;
countMatSafe(countMatSafe == 0) = 1;
S = S ./ repmat(countMatSafe, 1, 1, nf);
autoCountSafe = autoCount;
autoCountSafe(autoCountSafe == 0) = 1;
autoPower = autoPower ./ autoCountSafe;

if any(missMat(:))
    switch lower(params.method)
        case {'none', 'naive'}
            if params.verbose
                fprintf('[stitch_spectra_blocks] missing entries retained without completion\n');
            end
        case 'maxdet'
            if params.verbose
                fprintf('[stitch_spectra_blocks] running maxdet completion\n');
            end
            S = maxdet_completion_complex(S, missMat);
        case 'nnm'
            if params.verbose
                fprintf('[stitch_spectra_blocks] running NNM completion\n');
            end
            S = nnm_completion_complex(S, missMat);
        otherwise
            error('Unknown completion method: %s', params.method);
    end
end

if params.verbose
    if isscalar(params.lambda)
        lambdaText = sprintf('%.4g', params.lambda);
    else
        lambdaText = sprintf('[%.4g %.4g]', min(params.lambda(:)), max(params.lambda(:)));
    end
    fprintf('[stitch_spectra_blocks] regularizing inverse spectrum using %s (lambda=%s)\n', ...
        params.regularizer, lambdaText);
end
[P, Sreg] = glasso_precision_estimate(S, params.lambda, params);

meta.count_mat = rawCountMat;
meta.count_mat_safe = countMatSafe;
meta.missing_mask = missMat;
meta.normalize = params.normalize;
meta.auto_power = autoPower;
end


function Su = local_normalize_block(Su, normalizeMode)
mode = lower(string(normalizeMode));
switch mode
    case "none"
        return;
    case {"coherence", "coh"}
        n = size(Su, 1);
        nf = size(Su, 3);
        for f = 1:nf
            A = (Su(:,:,f) + Su(:,:,f)') / 2;
            p = real(diag(A));
            p(~isfinite(p) | p <= 0) = nan;
            denom = sqrt(p * p.');
            C = A ./ denom;
            C(~isfinite(C)) = 0;
            C(1:n+1:end) = 1;
            Su(:,:,f) = (C + C') / 2;
        end
    otherwise
        error('Unknown stitch normalize mode: %s', normalizeMode);
end
end
