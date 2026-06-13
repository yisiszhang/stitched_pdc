function resultOut = ibl_apply_glasso_to_result(resultIn, lambda, cfg, outFile)
%IBL_APPLY_GLASSO_TO_RESULT Apply original graphicalLasso to an existing stitched S.
%
% This bypasses re-stitching and uses resultIn.S directly. It is intended for
% fast glasso sensitivity runs after a baseline stitched result already exists.

if nargin < 4
    outFile = '';
end

S = resultIn.S;
freqs = resultIn.freqs(:);
[K, ~, nf] = size(S);
lambda = local_expand_lambda(lambda, nf);

if ~isfield(cfg, 'stitch')
    stitch = struct();
else
    stitch = cfg.stitch;
end
if ~isfield(stitch, 'glasso_maxiter') || isempty(stitch.glasso_maxiter)
    stitch.glasso_maxiter = 100;
end
if ~isfield(stitch, 'glasso_tol') || isempty(stitch.glasso_tol)
    stitch.glasso_tol = 1e-4;
end
if ~isfield(stitch, 'parallel') || isempty(stitch.parallel)
    stitch.parallel = false;
end
if ~isfield(stitch, 'verbose') || isempty(stitch.verbose)
    stitch.verbose = isfield(cfg, 'verbose') && cfg.verbose;
end

Sreg = zeros(size(S), 'like', S);
P = zeros(size(S), 'like', S);

if stitch.verbose
    fprintf('[ibl_apply_glasso_to_result] freqs=%d  lambda=[%.4g %.4g]  parallel=%d\n', ...
        nf, min(lambda), max(lambda), stitch.parallel);
end

if stitch.parallel
    parfor f = 1:nf
        [P(:,:,f), Sreg(:,:,f)] = local_glasso_slice(S(:,:,f), lambda(f), K, stitch);
    end
else
    for f = 1:nf
        if stitch.verbose
            fprintf('[ibl_apply_glasso_to_result] %4d/%4d  lambda=%.4g\n', f, nf, lambda(f));
        end
        [P(:,:,f), Sreg(:,:,f)] = local_glasso_slice(S(:,:,f), lambda(f), K, stitch);
    end
end

resultOut = resultIn;
resultOut.Sreg = Sreg;
resultOut.P = P;
resultOut.PDC = nonparam_pdc_H(Sreg, freqs, ...
    'maxiter', cfg.chronux.maxiter, 'tol', cfg.chronux.tol);
resultOut.glasso_lambda = lambda;
resultOut.glasso_source = 'ibl_apply_glasso_to_result';

if ~isempty(outFile)
    result = resultOut; %#ok<NASGU>
    save(outFile, 'result', '-v7.3');
end
end


function [Pfreq, SfreqReg] = local_glasso_slice(Sfreq, lambda, K, stitch)
S0 = (Sfreq + Sfreq') / 2;
Sreal = [real(S0) -imag(S0); imag(S0) real(S0)];
[~, W] = graphicalLasso(Sreal, lambda, stitch.glasso_maxiter, stitch.glasso_tol);
SfreqReg = W(1:K,1:K) + 1i * W(K+1:2*K,1:K);
SfreqReg = (SfreqReg + SfreqReg') / 2;
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
