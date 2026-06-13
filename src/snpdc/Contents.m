% SNPDC core functions.
%
% Main estimator
%   stitched_pdc              - Stitched non-parametric PDC from spectral blocks.
%   stitch_spectra_blocks     - Stitch/complete/regularize spectral blocks.
%   reconstruct_inversepsd    - Estimate spectra from time series and stitch.
%
% PDC and information flow
%   nonparam_pdc_H            - Non-parametric PDC via Wilson factorization.
%   pdc_to_infoflow           - Integrated information flow from PDC.
%
% Matrix completion and regularization
%   nnm_completion_complex    - Nuclear-norm completion for complex spectra.
%   maxdet_completion_complex - Max-det completion wrapper.
%   glasso_precision_estimate - Eigfloor/GLASSO inverse-spectrum regularizer.
%   graphicalLasso            - Graphical lasso solver.
%
% Spectral PCA
%   spectral_pca              - Frequency-wise PCA of Hermitian spectra.
%   spectral_pc_loading       - Fixed loading from band-averaged spectrum.
%   project_cross_spectrum    - Project cross spectra using fixed loadings.
%   point_spectral_pca        - Chronux point-process spectral PCA wrapper.
%
% IBL application helpers
%   ibl_default_config
%   ibl_scan_sessions
%   ibl_compute_area_pca
%   ibl_compute_cross_area_csd
%   ibl_stitch_saved_spectra
%   ibl_compare_allen_tracing
