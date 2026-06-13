function [vfix, dfix, Sbar] = ibl_fixed_spectral_loading(Saa, freqs, band)
%IBL_FIXED_SPECTRAL_LOADING Backward-compatible alias for SPECTRAL_PC_LOADING.
[vfix, dfix, Sbar] = spectral_pc_loading(Saa, freqs, band);
end
