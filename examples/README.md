# abqhom examples

Each example is a self-contained `main.py`. Run it from its own folder with

```
abaqus cae noGUI=main.py
```

(e.g. `abq2024le cae noGUI=main.py` for the Learning Edition). The abqhom
source is located automatically relative to the example; only the `gmsh`
(and optional `Paraqus`, for the VTK export) path needs to be set at the
top of each script. Results are written next to the script into `abq/`
(Abaqus job files), `results/` (the effective stiffness as CSV) and
`vtk_output/` (Paraqus VTK files).

## Beam / truss lattices

| Example | Cell | Check |
|---|---|---|
| `lattice_simple_cubic` | edge struts of the cube | analytical `C11 = E*pi*r^2/L^2` |
| `lattice_bcc` | cubic edges + body diagonals | cubic symmetry |
| `lattice_honeycomb_2d` | regular hexagonal honeycomb | in-plane isotropy |
| `lattice_kelvin_3d` | Kelvin (truncated octahedron) | cubic symmetry |
| `lattice_plasticity` | simple cubic, elasto-plastic | nonlinear (explicit) |

The lattices whose struts lie on the bounding box (simple cubic, BCC,
Kelvin) rely on the automatic periodic sharing weights in
`average_stress`; no manual stiffness correction is applied.

## Continuum and mixed RVEs

| Example | RVE | Check |
|---|---|---|
| `fiber_reinforced_2d` | fibers in a matrix (plane strain) | Voigt / Reuss bounds |
| `spherical_inclusion_3d` | spherical inclusion (manual boundary) | Voigt / Reuss bounds |
| `beam_fibers_in_matrix_3d` | beam fibers embedded in a continuum | transverse isotropy |
