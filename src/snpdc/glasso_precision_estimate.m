function [P, Sreg] = glasso_precision_estimate(S, lambda, opts)
%GLASSO_PRECISION_ESTIMATE Regularize inverse spectrum for downstream Wilson/PDC.
%
%   lambda can be scalar or an nf-by-1 vector of frequency-specific penalties.

if nargin < 3 || isempty(opts)
    opts = struct();
end
if ~isfield(opts, 'regularizer') || isempty(opts.regularizer)
    opts.regularizer = 'eigfloor';
end
if ~isfield(opts, 'min_eig_rel') || isempty(opts.min_eig_rel)
    opts.min_eig_rel = 1e-6;
end
if ~isfield(opts, 'min_eig_abs') || isempty(opts.min_eig_abs)
    opts.min_eig_abs = 1e-8;
end
if ~isfield(opts, 'glasso_fallback_lambda') || isempty(opts.glasso_fallback_lambda)
    opts.glasso_fallback_lambda = 0.01;
end
if ~isfield(opts, 'glasso_maxiter') || isempty(opts.glasso_maxiter)
    opts.glasso_maxiter = 30;
end
if ~isfield(opts, 'glasso_tol') || isempty(opts.glasso_tol)
    opts.glasso_tol = 1e-4;
end
if ~isfield(opts, 'glasso_retry_factor') || isempty(opts.glasso_retry_factor)
    opts.glasso_retry_factor = 5;
end
if ~isfield(opts, 'glasso_standardize') || isempty(opts.glasso_standardize)
    opts.glasso_standardize = true;
end
if ~isfield(opts, 'parallel') || isempty(opts.parallel)
    opts.parallel = false;
end
if ~isfield(opts, 'glasso_diagnostics') || isempty(opts.glasso_diagnostics)
    opts.glasso_diagnostics = false;
end
if ~isfield(opts, 'glasso_live_print') || isempty(opts.glasso_live_print)
    opts.glasso_live_print = false;
end
if ~isfield(opts, 'verbose') || isempty(opts.verbose)
    opts.verbose = false;
end

[~, K, nf] = size(S);
lambda = local_expand_lambda(lambda, nf);

tolEig = 1e-12;
regMode = lower(string(opts.regularizer));

Sreg = zeros(size(S), 'like', S);
P = zeros(size(S), 'like', S);

needsRepairVec = false(nf, 1);
usedGlassoVec = false(nf, 1);
usedFallbackVec = false(nf, 1);
usedRetryVec = false(nf, 1);
lambdaUsed = nan(nf, 1);
hitMaxIterVec = false(nf, 1);
hitShootingMaxIterVec = false(nf, 1);
glassoIterVec = nan(nf, 1);

if opts.parallel
    parOpts = opts;
    parOpts.verbose = false;
    parfor i = 1:nf
        [P(:,:,i), Sreg(:,:,i), needsRepairVec(i), usedGlassoVec(i), ...
            usedFallbackVec(i), usedRetryVec(i), lambdaUsed(i), ...
            hitMaxIterVec(i), hitShootingMaxIterVec(i), glassoIterVec(i)] = ...
            local_process_frequency(S(:,:,i), K, lambda(i), regMode, tolEig, parOpts, i);
    end
else
    for i = 1:nf
        [P(:,:,i), Sreg(:,:,i), needsRepairVec(i), usedGlassoVec(i), ...
            usedFallbackVec(i), usedRetryVec(i), lambdaUsed(i), ...
            hitMaxIterVec(i), hitShootingMaxIterVec(i), glassoIterVec(i)] = ...
            local_process_frequency(S(:,:,i), K, lambda(i), regMode, tolEig, opts, i);
        if opts.glasso_live_print && (regMode == "glasso" || regMode == "auto")
            local_print_frequency_status(i, nf, lambdaUsed(i), usedGlassoVec(i), ...
                usedRetryVec(i), usedFallbackVec(i), hitMaxIterVec(i), ...
                hitShootingMaxIterVec(i), glassoIterVec(i));
        end
    end
end

nBad = sum(needsRepairVec & ~usedGlassoVec);
nGlasso = sum(usedGlassoVec);
nFallback = sum(usedFallbackVec);
nRetry = sum(usedRetryVec);

if opts.verbose
    switch regMode
        case "glasso"
            fprintf('[glasso_precision_estimate] glasso slices=%d  retry=%d  eigfloor-fallback=%d\n', ...
                nGlasso, nRetry, nFallback);
            fprintf('[glasso_precision_estimate] lambda used range=[%.4g %.4g]\n', ...
                min(lambdaUsed), max(lambdaUsed));
            if opts.glasso_diagnostics
                local_print_maxiter_diagnostics(hitMaxIterVec, hitShootingMaxIterVec, glassoIterVec);
            end
        case "auto"
            fprintf('[glasso_precision_estimate] eigfloor-repairs=%d  glasso-slices=%d  retry=%d  eigfloor-fallback=%d\n', ...
                nBad, nGlasso, nRetry, nFallback);
        otherwise
            fprintf('[glasso_precision_estimate] eigfloor repairs=%d/%d slices\n', nBad, nf);
    end
end
end


function [Pfreq, SfreqReg, needsRepair, usedGlasso, usedFallback, usedRetry, usedLambda, ...
    hitMaxIter, hitShootingMaxIter, glassoIter] = ...
    local_process_frequency(Sfreq, K, lam, regMode, tolEig, opts, freqIdx)
A = local_hermitian(Sfreq);
d = eig(A, 'vector');
needsRepair = ~all(isfinite(d)) || min(real(d)) < tolEig;
hitMaxIter = false;
hitShootingMaxIter = false;
glassoIter = nan;

switch regMode
    case "glasso"
        [Areg, usedGlasso, usedFallback, usedRetry, usedLambda, hitMaxIter, hitShootingMaxIter, glassoIter] = ...
            local_glasso_or_fallback(A, K, lam, opts, freqIdx);
    case "auto"
        if lam > 0 && ~needsRepair
            [Areg, usedGlasso, usedFallback, usedRetry, usedLambda, hitMaxIter, hitShootingMaxIter, glassoIter] = ...
                local_glasso_or_fallback(A, K, lam, opts, freqIdx);
        else
            Areg = local_project_spd(A, opts.min_eig_rel, opts.min_eig_abs);
            usedGlasso = false;
            usedFallback = false;
            usedRetry = false;
            usedLambda = lam;
        end
    otherwise
        if needsRepair
            Areg = local_project_spd(A, opts.min_eig_rel, opts.min_eig_abs);
        else
            Areg = A;
        end
        usedGlasso = false;
        usedFallback = false;
        usedRetry = false;
        usedLambda = lam;
end

SfreqReg = local_hermitian(Areg);
Pfreq = SfreqReg \ eye(K);
end


function lambda = local_expand_lambda(lambda, nf)
if isscalar(lambda)
    lambda = repmat(lambda, nf, 1);
else
    lambda = lambda(:);
    assert(numel(lambda) == nf, ...
        'Frequency-specific lambda must have one value per frequency bin.');
end
lambda(~isfinite(lambda) | lambda < 0) = 0;
end


function [Areg, usedGlasso, usedFallback, usedRetry, usedLambda, hitMaxIter, hitShootingMaxIter, glassoIter] = ...
    local_glasso_or_fallback(A, K, lambda, opts, freqIdx)
usedGlasso = false;
usedFallback = false;
usedRetry = false;
usedLambda = lambda;
hitMaxIter = false;
hitShootingMaxIter = false;
glassoIter = nan;

if lambda <= 0
    Areg = local_project_spd(A, opts.min_eig_rel, opts.min_eig_abs);
    usedFallback = true;
    return;
end

Sreal = [real(A) -imag(A); imag(A) real(A)];
[Swork, scaleVec] = local_standardize_real_cov(Sreal, opts.glasso_standardize);

[ok, W, info] = local_try_glasso(Swork, lambda, opts);
hitMaxIter = hitMaxIter || local_info_flag(info, 'hit_maxiter');
hitShootingMaxIter = hitShootingMaxIter || local_info_flag(info, 'n_shooting_maxed');
glassoIter = local_info_value(info, 'iterations');
if ~ok
    retryLambda = max([lambda * opts.glasso_retry_factor, opts.glasso_fallback_lambda, lambda]);
    if isfinite(retryLambda) && retryLambda > lambda
        usedRetry = true;
        if opts.verbose
            fprintf('[glasso] freq=%d retry lambda %.4g -> %.4g\n', freqIdx, lambda, retryLambda);
        end
        [ok, W, info] = local_try_glasso(Swork, retryLambda, opts);
        hitMaxIter = hitMaxIter || local_info_flag(info, 'hit_maxiter');
        hitShootingMaxIter = hitShootingMaxIter || local_info_flag(info, 'n_shooting_maxed');
        glassoIter = local_info_value(info, 'iterations');
        usedLambda = retryLambda;
    end
end

if ok
    W = local_unstandardize_real_cov(W, scaleVec);
    Areg = local_complex_from_real_cov(W, K);
    Areg = local_project_spd(Areg, opts.min_eig_rel, opts.min_eig_abs);
    usedGlasso = true;
    return;
end

Areg = local_project_spd(A, opts.min_eig_rel, opts.min_eig_abs);
usedFallback = true;
if opts.verbose
    fprintf('[glasso] freq=%d fallback to eigfloor (lambda=%.4g)\n', freqIdx, usedLambda);
end
end


function [R, scaleVec] = local_standardize_real_cov(Sreal, doStandardize)
scaleVec = ones(size(Sreal,1), 1);
if ~doStandardize
    R = Sreal;
    return;
end
d = real(diag(Sreal));
valid = isfinite(d) & d > 0;
if any(valid)
    d(~valid) = median(d(valid));
else
    d(:) = 1;
end
scaleVec = sqrt(d(:));
R = Sreal ./ (scaleVec * scaleVec.');
R = (R + R') / 2;
R(1:size(R,1)+1:end) = 1;
end


function Sreal = local_unstandardize_real_cov(R, scaleVec)
Sreal = R .* (scaleVec * scaleVec.');
Sreal = (Sreal + Sreal') / 2;
end


function [ok, W, info] = local_try_glasso(Sreal, rho, opts)
ok = false;
W = [];
info = struct();
try
    if opts.glasso_diagnostics
        [~, W, info] = graphicalLassoDiag(Sreal, rho, opts.glasso_maxiter, opts.glasso_tol);
    else
        [~, W] = graphicalLasso(Sreal, rho, opts.glasso_maxiter, opts.glasso_tol);
    end
catch
    return;
end
ok = all(isfinite(W(:)));
end


function tf = local_info_flag(info, field)
tf = false;
if isempty(info) || ~isfield(info, field)
    return;
end
val = info.(field);
if islogical(val)
    tf = any(val(:));
else
    tf = any(val(:) ~= 0);
end
end


function val = local_info_value(info, field)
val = nan;
if ~isempty(info) && isfield(info, field)
    val = info.(field);
end
end


function local_print_maxiter_diagnostics(hitOuter, hitInner, iterVec)
outerIdx = find(hitOuter);
innerIdx = find(hitInner);
fprintf('[glasso_precision_estimate] outer maxIt frequencies: %s\n', local_idx_text(outerIdx));
fprintf('[glasso_precision_estimate] shooting maxIt frequencies: %s\n', local_idx_text(innerIdx));
if any(isfinite(iterVec))
    fprintf('[glasso_precision_estimate] glasso iteration range=[%g %g]\n', ...
        min(iterVec(isfinite(iterVec))), max(iterVec(isfinite(iterVec))));
end
end


function local_print_frequency_status(i, nf, lambdaUsed, usedGlasso, usedRetry, usedFallback, hitOuter, hitInner, iterCount)
if isnan(iterCount)
    iterText = 'n/a';
else
    iterText = sprintf('%g', iterCount);
end
fprintf('[glasso freq] %4d/%4d  lambda=%.4g  glasso=%d  retry=%d  fallback=%d  outerMax=%d  shootingMax=%d  iter=%s\n', ...
    i, nf, lambdaUsed, usedGlasso, usedRetry, usedFallback, hitOuter, hitInner, iterText);
end


function txt = local_idx_text(idx)
if isempty(idx)
    txt = 'none';
elseif numel(idx) <= 30
    txt = mat2str(idx(:)');
else
    txt = sprintf('%d frequencies, first 30=%s', numel(idx), mat2str(idx(1:30)'));
end
end


function A = local_complex_from_real_cov(W, K)
A = W(1:K,1:K) + 1i * W(K+1:2*K,1:K);
A = local_hermitian(A);
end


function Areg = local_project_spd(A, minEigRel, minEigAbs)
A = local_hermitian(A);
[V, D] = eig(A, 'vector');
d = real(D);

scale = real(trace(A)) / max(size(A,1), 1);
if ~isfinite(scale) || scale <= 0
    scale = 1;
end
floorVal = max(minEigAbs, minEigRel * scale);
d(~isfinite(d)) = floorVal;
d(d < floorVal) = floorVal;

Areg = V * diag(d) * V';
Areg = local_hermitian(Areg);
end


function A = local_hermitian(A)
A = (A + A') / 2;
end
