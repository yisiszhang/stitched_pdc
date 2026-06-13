function sab = ibl_project_cross_spectrum(Sab, va, vb)
%IBL_PROJECT_CROSS_SPECTRUM Project neuron-level cross-spectrum to area PC1s.

va = squeeze(va);
vb = squeeze(vb);
nFreq = size(Sab, 3);
sab = zeros(1, nFreq);

for f = 1:nFreq
    sab(f) = va(:,f)' * Sab(:,:,f) * vb(:,f);
end
end
