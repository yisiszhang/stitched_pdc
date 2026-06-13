function [evals, evecs, explVar] = spectral_pca(Saa, nPcs)
%SPECTRAL_PCA Frequency-wise PCA of a Hermitian cross-spectrum.
%
%   [evals, evecs, explVar] = spectral_pca(Saa, nPcs)
%
% Inputs
%   Saa  : nSignals-by-nSignals-by-nFreq auto/cross-spectrum within one area
%   nPcs : number of PCs to retain
%
% Outputs
%   evals   : nPcs-by-nFreq eigenvalues
%   evecs   : nSignals-by-nPcs-by-nFreq PC loadings
%   explVar : nPcs-by-nFreq explained variance fractions
%
% Sign convention
%   Each loading is phase-aligned so its largest-magnitude entry is real
%   positive. This removes arbitrary eigensolver phase/sign flips.

[n, ~, nFreq] = size(Saa);
nPcs = min(nPcs, n);

evals = zeros(nPcs, nFreq);
evecs = zeros(n, nPcs, nFreq);

for f = 1:nFreq
    Sf = (Saa(:,:,f) + Saa(:,:,f)') / 2;
    Sf(~isfinite(Sf)) = 0;
    [V, D] = eig(Sf, 'vector');
    [d, ord] = sort(real(D), 'descend');
    d = max(d, 0);
    V = V(:, ord);
    for pc = 1:nPcs
        V(:, pc) = orient_loading(V(:, pc));
    end
    evals(:,f) = d(1:nPcs);
    evecs(:,:,f) = V(:, 1:nPcs);
end

traceVals = squeeze(real(sum(local_diag3(Saa), 1)));
traceVals(~isfinite(traceVals)) = NaN;
traceVals(traceVals <= 0) = 1;
explVar = evals ./ traceVals;
end


function d = local_diag3(S)
n = size(S,1);
nFreq = size(S,3);
d = zeros(n, nFreq);
for k = 1:nFreq
    d(:,k) = diag(S(:,:,k));
end
end
