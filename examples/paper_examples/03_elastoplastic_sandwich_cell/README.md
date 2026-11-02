# Example 3: elastoplastic sandwich cell

This example reproduces the cyclic and yield-surface studies reported in
Section 5.3 of the accompanying manuscript. The mixed-dimensional cell
contains a soft elastic continuum matrix, two ideally plastic shell skins,
and an embedded ideally plastic BCC beam core.

## Scripts

- `main.py` computes the cyclic effective response, PEEQ histories,
  verification measures, and VTK output.
- `yield_surface_sections.py` computes the elastic calibration and the
  initial and evolved yield-surface sections in the `sigma11-sigma22` and
  `sigma11-sigma33` planes.
- `yield_surface_3d.py` samples the complete initial or evolved yield surface
  in normal-stress space using independent resumable workers.
- `sandwich.py` provides the shared Gmsh geometry definition.

The scripts are compatible with both Python 2.7-based legacy Abaqus releases
and current Python 3-based releases. Activate the isolated Abaqus runtime as
described in the project documentation.

## Cyclic response

Run the cyclic analysis directly from this directory:

```text
abaqus cae noGUI=main.py
```

The prescribed axial history contains two fully reversed cycles between
`epsilon_11 = -0.03` and `epsilon_11 = 0.03`, followed by unloading to zero.

## Yield-surface sections

Run the planar sampling with:

```text
abaqus cae noGUI=yield_surface_sections.py
```

The default configuration computes the elastic calibration and samples both
planes for the initial and evolved states with 36 directions per plane. The
configuration block accepts `selected_state` equal to `initial`, `evolved`,
or `both`, and `selected_plane` equal to `sigma11_sigma22`,
`sigma11_sigma33`, or `None`. Set `calibration_only=True` when only the
effective stiffness and compliance needed by the 3D sampler are required.
Completed angular prefixes are resumed from their existing CSV files.

## Three-dimensional yield surfaces

The 3D sampler requires the elastic calibration produced by
`yield_surface_sections.py`. Configure `state_name`, `worker_index`, and
`number_of_workers` in `yield_surface_3d.py`, then run every worker with:

```text
abaqus cae noGUI=yield_surface_3d.py
```

Worker indices start at zero and must cover
`0, ..., number_of_workers - 1`. All workers for one state must use the same
worker count and analysis configuration. Each worker resumes its assigned
subset of the 642 sampling directions independently. After all workers for a
state have completed, set `merge_only=True` and run the script once more with
the same `state_name` and `number_of_workers`. Repeat the worker and merge
procedure for the other state when both surfaces are required.

Configuration records accompany the resumable yield-surface results. If a
material, mesh, loading, or sampling setting changes, remove the corresponding
generated result directory before restarting that analysis.

## Generated output

The cyclic analysis writes:

| Directory | Contents |
|---|---|
| `abq/` | Cyclic Abaqus job |
| `results/` | Effective response, PEEQ history, and verification measures |
| `vtk_output/` | Cyclic VTK collection |

The cyclic `results/` directory contains:

| File | Contents |
|---|---|
| `cyclic_response.csv` | Pseudo-time, axial strain, nine effective-stress components, stress-asymmetry residual, and Hill–Mandel residual |
| `peeq_history.csv` | Pseudo-time and maximum shell and beam PEEQ |
| `verification_metrics.csv` | Maximum stress asymmetry, Hill–Mandel residual, shell PEEQ, and beam PEEQ |

Yield-surface output is stored below `yield_surface/`:

| Directory | Contents |
|---|---|
| `calibration_abq/` | Elastic-calibration Abaqus jobs |
| `results/` | Elastic calibration and initial yield-surface sections |
| `evolved_results/` | Evolved yield-surface sections |
| `initial_abq/`, `evolved_abq/` | Section-sampling Abaqus jobs |
| `surface_3d/initial/`, `surface_3d/evolved/` | Worker files and merged 3D surfaces |

The sectional result directories contain `sigma11_sigma22.csv`,
`sigma11_sigma33.csv`, the common base strain and stress, and their sampling
configuration. Each merged 3D directory contains `yield_surface.csv`, the
sampling directions in `directions.csv`, triangular connectivity in
`faces.csv`, the common base state, and the resumable worker files below
`parts/`.
