# abqhom examples

The examples are divided into concise user tutorials and the more extensive
manuscript-reproduction cases in `paper_examples`. Each tutorial is a
self-contained `main.py`. Run it from its own folder with

```
abaqus cae noGUI=main.py
```

The abqhom source is located automatically relative to the example. Before
launching Abaqus, add the directories containing the Gmsh and optional
Paraqus Python modules to `PYTHONPATH`. Results are written next to the script
into `abq/` (Abaqus job files), `results/` (CSV data and plots) and, when
Paraqus is available, `vtk_output/` (VTK files). VTK export is skipped when
Paraqus is not installed.

For example, in PowerShell:

```powershell
$env:PYTHONPATH = (
    "C:/path/to/gmsh/lib;" +
    "C:/path/to/Paraqus/src;" +
    $env:PYTHONPATH
)
```

## Beam / truss lattices

| Example | Cell | Purpose |
|---|---|---|
| `lattice_simple_cubic` | edge struts of the cube | analytical `C11 = E*pi*r^2/L^2` |
| `lattice_honeycomb_2d` | regular hexagonal honeycomb | in-plane isotropy |
| `lattice_kelvin_3d` | Kelvin (truncated octahedron) | finite-strain beam response |
| `lattice_plasticity` | simple cubic, cyclic perfect plasticity | hysteresis |

The simple-cubic and Kelvin cells rely on the automatic removal of redundant
secondary element copies. No manual stiffness correction is applied.

## Continuum and mixed RVEs

| Example | RVE | Purpose |
|---|---|---|
| `fiber_reinforced_2d` | fibers in a matrix (plane strain) | Voigt / Reuss bounds |
| `spherical_inclusion_3d` | spherical inclusion (manual boundary) | finite-strain two-phase response |
| `beam_fibers_in_matrix_3d` | beam fibers embedded in a continuum | mixed-dimensional embedded coupling |

## Shell RVEs

| Example | RVE | Purpose |
|---|---|---|
| `auxetic_lattice` | re-entrant honeycomb, shell walls | negative in-plane Poisson's ratio |
| `gyroid_shell_3d` | sheet-gyroid TPMS with an imported shell mesh | finite-strain shell response |

The `auxetic_lattice` cell is the classic 2D re-entrant (bow-tie) honeycomb
given a finite depth: its ribs (in the x-z plane) are extruded along y into
thin plate walls. The inclined ribs point inward (re-entrant angle), so
stretching the honeycomb in-plane swings them outward and the cell expands
laterally -- a negative Poisson's ratio in the x-z plane. Along the extrusion
axis y the walls just stretch, so the effective stiffness is orthotropic
(auxetic in-plane, an ordinary stiff honeycomb prism along y). The in-plane
honeycomb is built by clipping the infinite lattice to one period before
extrusion, like the beam honeycomb example.

The gyroid is an implicit surface. Its periodic surface mesh is generated with
the [microgen](https://github.com/3MAH/microgen) package: it is pre-computed
and exported to `gyroid.msh` by `generate_mesh.py` (run once in an
environment that has microgen installed), so Abaqus itself does not need
microgen -- `main.py` only reads the mesh file and hands it to abqhom via
`GmshRVE.read_from_mesh`.

## Manuscript-reproduction cases

The three cases under `paper_examples` reproduce the numerical examples of
the accompanying manuscript. They include multi-stage analyses and plotting
scripts and are therefore intentionally more detailed than the tutorials.
See `paper_examples/README.md` for their execution order and purpose.
