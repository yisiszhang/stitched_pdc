# Stitched Non-Parametric PDC

MATLAB code for estimating non-parametric partial directed coherence (PDC)
from partially observed, overlapping recordings. The package stitches local
cross-spectral density blocks into a global spectrum, optionally normalizes
blocks to coherence, completes missing entries, regularizes the spectrum, and
computes PDC/information flow.

Direction convention:

- `PDC(i,j,f)` is influence from source `j` to target `i`.
- `infoflow(i,j)` is integrated information flow from source `j` to target `i`.

## Repository Layout

- `src/snpdc/`: core algorithms and reusable data-analysis helpers.
- `examples/`: lightweight runnable demos and data-application templates.
- `external/`: bundled third-party code; see `THIRD_PARTY_NOTICES.md`.
- `docs/`: API and workflow notes.
- `tests/`: smoke-test area for future tests.

Generated outputs such as `ibl_output/` are intentionally not part of the
public package interface.

## Quickstart

Open MATLAB in this folder:

```matlab
setup
run(fullfile('examples', 'var3node_demo.m'))
```

Optional MatrixCompletion MEX compilation:

```matlab
setup compile
```

The code runs without the optional MEX, but nuclear-norm matrix completion can
be slower.

## Minimal Core Usage

If you already have local spectral matrices:

```matlab
setup

% Observed blocks: Sblocks{u} is |recset{u}| x |recset{u}| x nFreq.
recset = {[1 2 3], [3 4 5], [5 6 7]};
f = f(:);

params.method = 'nnm';          % matrix completion
params.normalize = 'coherence'; % useful across heterogeneous recordings
params.regularizer = 'eigfloor';
params.lambda = 0;
params.pdc_maxiter = 100;
params.pdc_tol = 1e-3;

result = stitched_pdc(Sblocks, recset, f, params);
imagesc(result.infoflow);
xlabel('Source');
ylabel('Target');
```

## Spectral PCA for Multi-Neuron Areas

Reusable spectral-PCA functions are included for point-process/area analyses:

- `spectral_pca`: frequency-wise PCA of a Hermitian spectrum.
- `spectral_pc_loading`: fixed PC loading from a band-averaged spectrum.
- `project_cross_spectrum`: project a multi-neuron cross-spectrum to one
  area-level cross-spectrum.
- `point_spectral_pca`: Chronux point-process wrapper for spike times.

The PC sign/phase convention is fixed by rotating each loading so the
largest-magnitude element is real positive.


It requires:

- Chronux with `CrossSpecMatpt` on the MATLAB path.
- Editing only the configuration block at the top of the script.

See `docs/WORKFLOWS.md` for the step-by-step pipeline.

## Documentation

- `docs/API_REFERENCE.md`
- `docs/WORKFLOWS.md`
- `examples/README.md`

## Dependencies

- MATLAB R2020b+ recommended.
- Signal Processing Toolbox for `cpsd`-based continuous spectra.
- Chronux for point-process multitaper spectra and the IBL workflow.
- Parallel Computing Toolbox is optional.

## Third-Party Code

Third-party code is included in `external/` for convenience. See
`THIRD_PARTY_NOTICES.md` for licenses and provenance.
