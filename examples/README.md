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
| `lattice_plasticity` | simple cubic, cyclic elasto-plastic | isotropic hardening, hysteresis |

The lattices whose struts lie on the bounding box (simple cubic, BCC,
Kelvin) rely on the automatic periodic sharing weights in
`average_stress`; no manual stiffness correction is applied.

## Continuum and mixed RVEs

| Example | RVE | Check |
|---|---|---|
| `fiber_reinforced_2d` | fibers in a matrix (plane strain) | Voigt / Reuss bounds |
| `spherical_inclusion_3d` | spherical inclusion (manual boundary) | Voigt / Reuss bounds |
| `beam_fibers_in_matrix_3d` | beam fibers embedded in a continuum | transverse isotropy |

## Shell RVEs

| Example | RVE | Check |
|---|---|---|
| `auxetic_lattice` | re-entrant honeycomb, shell walls | negative in-plane Poisson's ratio |
| `gyroid_shell_microgen_3d` | sheet-gyroid TPMS (shell), meshed with microgen | cubic symmetry |

The `auxetic_lattice` cell is the classic 2d re-entrant (bow-tie) honeycomb
given a finite depth: its ribs (in the x-z plane) are extruded along y into
thin plate walls. The inclined ribs point inward (re-entrant angle), so
stretching the honeycomb in-plane swings them outward and the cell expands
laterally -- a negative Poisson's ratio in the x-z plane. Along the extrusion
axis y the walls just stretch, so the effective stiffness is orthotropic
(auxetic in-plane, an ordinary stiff honeycomb prism along y). The in-plane
honeycomb is built by clipping the infinite lattice to one period before
extrusion, like the beam honeycomb example.

The gyroid is an implicit surface, which Gmsh cannot build as CAD
geometry. The periodic surface mesh is generated with the
[microgen](https://github.com/3MAH/microgen) package: it is pre-computed
and exported to `gyroid.msh` by `generate_mesh.py` (run once in an
environment that has microgen installed), so Abaqus itself does not need
microgen -- `main.py` only reads the mesh file and hands it to abqhom via
`GmshRVE.read_from_mesh`.
