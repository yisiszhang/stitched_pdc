# Examples

Run `setup` from the repository root before running examples.

## Self-contained demos

- `var3node_demo.m`: three-node VAR example showing how observing the full
  network mitigates spurious pairwise PDC.
- `izh3node_pairwise_demo.m`: three Izhikevich neurons with stitched pair
  observations and information-flow comparison.
- `calcium_imaging_demo.m`: synthetic calcium-imaging-like traces with
  overlapping fields of view.

## Data application templates

- `ibl_spontaneous_pipeline.m`: compact IBL Neuropixels spontaneous-activity
  pipeline template. Requires a local OpenAlyx/ONE cache and Chronux.
- `ibl_within_session_validation_demo.m`: within-session subset-stitching
  validation template for IBL data.

Large manuscript reproductions are in `scripts/`, not `examples/`, so the
example folder remains runnable and lightweight.
