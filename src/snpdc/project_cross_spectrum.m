function sab = project_cross_spectrum(Sab, va, vb)
%PROJECT_CROSS_SPECTRUM Project a cross-spectrum with two fixed loadings.
%
%   sab = project_cross_spectrum(Sab, va, vb)
%
% Inputs
%   Sab : nA-by-nB-by-nFreq cross-spectrum between two multineuron areas
%   va  : nA-by-1 loading for area A
%   vb  : nB-by-1 loading for area B
%
% Output
%   sab : 1-by-nFreq projected area-level cross-spectrum

va = va(:);
vb = vb(:);
nFreq = size(Sab, 3);
sab = zeros(1, nFreq);

for f = 1:nFreq
    sab(f) = va' * Sab(:,:,f) * vb;
end
end
