function result = stitched_pdc(Sblocks, recset, f, params)
%STITCHED_PDC Core stitched non-parametric PDC estimator.
%
%   result = stitched_pdc(Sblocks, recset, f, params)
%
% Inputs
%   Sblocks : cell array of observed block spectra. Sblocks{u} is
%             |Omega_u|-by-|Omega_u|-by-nFreq.
%   recset  : cell array of global node indices for each observed block.
%   f       : nFreq-by-1 frequency grid.
%   params  : struct controlling completion, normalization, regularization,
%             and Wilson factorization.
%
% Common params fields
%   method      : 'none'|'naive'|'nnm'|'maxdet' (default: 'nnm')
%   normalize   : 'none'|'coherence' (default: 'none')
%   regularizer : 'eigfloor'|'glasso' (default: 'eigfloor')
%   lambda      : scalar or nFreq-vector GLASSO lambda (default: 0)
%   pdc_maxiter : Wilson iterations (default: 100)
%   pdc_tol     : Wilson tolerance (default: 1e-3)
%
% Output fields
%   PDC      : K-by-K-by-nFreq non-parametric PDC
%   infoflow : integrated information flow, column -> row
%   S        : stitched spectrum before regularization
%   Sreg     : spectrum after eigfloor/GLASSO regularization
%   P        : inverse spectrum
%   meta     : observation counts and missing-pair metadata

arguments
    Sblocks (1,:) cell
    recset (1,:) cell
    f (:,1) double
    params struct = struct()
end

if ~isfield(params, 'method'); params.method = 'nnm'; end
if ~isfield(params, 'normalize'); params.normalize = 'none'; end
if ~isfield(params, 'regularizer'); params.regularizer = 'eigfloor'; end
if ~isfield(params, 'lambda'); params.lambda = 0; end
if ~isfield(params, 'pdc_maxiter'); params.pdc_maxiter = 100; end
if ~isfield(params, 'pdc_tol'); params.pdc_tol = 1e-3; end

[P, S, Sreg, f, meta] = stitch_spectra_blocks(Sblocks, recset, f, params);
PDC = nonparam_pdc_H(Sreg, f, 'maxiter', params.pdc_maxiter, 'tol', params.pdc_tol);

result = struct();
result.PDC = PDC;
result.infoflow = pdc_to_infoflow(PDC);
result.S = S;
result.Sreg = Sreg;
result.P = P;
result.freqs = f;
result.recset = recset;
result.params = params;
result.meta = meta;
end
