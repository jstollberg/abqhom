# Elastoplastic sandwich cell

This example reproduces the cyclic and yield-surface analyses of the
elastoplastic sandwich cell. The continuum matrix is elastic, while the shell
skins and embedded beam core obey ideal J2 plasticity.

## Scripts

- `main.py` computes the cyclic effective response, verification measures,
  PEEQ history, and optional VTK output.
- `yield_surface_sections.py` computes the effective elastic calibration and
  the initial and evolved yield-surface sections.
- `yield_surface_3d.py` samples the complete initial or evolved
  three-dimensional yield surface using independent resumable workers and
  merges their results.
- `sandwich.py` contains the shared Gmsh geometry definition.

## Execution order

The cyclic analysis is independent and can be run directly with Abaqus. For
the yield surfaces, first run `yield_surface_sections.py`. Its elastic
calibration provides the compliance used to map stress-space directions to
strain-controlled probes. Then configure and run the required workers in
`yield_surface_3d.py`. After all workers have completed, set `merge_only` to
`True` and run the 3D script once more for each state.

All analysis settings are defined in the configuration blocks of the entry
scripts. Resumable outputs include configuration records. If a material,
mesh, loading, or sampling setting changes, remove the corresponding generated
results before restarting that analysis. Generated Abaqus files, CSV data, and
VTK output are ignored by Git.
