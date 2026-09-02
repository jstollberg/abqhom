# abqhom

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

A small toolbox for linear (and simple nonlinear) computational
homogenization of periodic representative volume elements (RVEs) with
[Gmsh](https://gmsh.info/) and [Abaqus](https://www.3ds.com/products/simulia/abaqus).
You build a unit cell in Gmsh, and `abqhom` prepares the periodic model,
runs the analysis in Abaqus and returns the effective (homogenized)
stiffness.

## Features

- **Beam, truss, shell and continuum RVEs** — and mixtures of them in one
  model. Beam and shell sections are integrated during the analysis, so
  nonlinear material behaviour is captured.
- **General unit cells** — the struts/walls need not connect to the
  corners of the bounding box. Classic "cubic" cells as well as
  honeycomb, Kelvin and other lattices are supported; the boundary and
  its periodic connections are detected automatically from the geometry
  (or supplied manually). A ready-made mesh (e.g. an implicit TPMS surface
  such as a gyroid, meshed outside the CAD kernel) can be ingested
  directly with `GmshRVE.read_from_mesh`, which classifies its boundary
  nodes and checks it for periodicity.
- **Periodic and affine (Dirichlet) displacement** boundary conditions,
  with automatic suppression of rigid-body motion. Periodic BCs support
  beam/shell rotational degrees of freedom; the affine BCs are for
  continua only (they raise for rotational/warping DOFs).
- **Effective stress from element nodal forces (`NFORC`)** via the
  discrete Hill–Mandel identity — one material-agnostic routine for
  periodic and displacement BCs, for continuum, beam/truss and shell
  elements.
- **Automatic periodic strut-sharing weights**: struts lying on a
  bounding-box face (shared by 2 cells) or edge (shared by 4) contribute
  only their in-cell fraction to the effective stiffness.
- **VTK export** of the results via
  [Paraqus](https://github.com/InstituteOfMechanics/Paraqus) (optional),
  with per-strut radii taken automatically from the assigned sections.

## Requirements

- **Abaqus** (CAE and the solver). Scripts that use `AbaqusRVE` must be
  run in the Abaqus Python interpreter, e.g. `abaqus cae noGUI=main.py`.
  Tested with Abaqus 2024 LE (Python 3.10) and Abaqus 2023 LE
  (Python 2.7); the code is Python 2/3 compatible.
- **Gmsh** with its Python API (the SDK from https://gmsh.info/).
- **NumPy** (shipped with Abaqus).
- **Paraqus** (optional, only for the VTK export).

## Installation

`abqhom` is used directly from source — put `src` on the Python path so
that `import abqhom` works, and make sure Gmsh is importable:

```python
import sys
sys.path.insert(0, "/path/to/abqhom/src")
sys.path.append("/path/to/gmsh/lib")        # folder containing gmsh.py
# optional, for the VTK export:
sys.path.append("/path/to/Paraqus/src")
```

The examples locate the `abqhom` source relative to their own location
automatically; only the Gmsh (and optional Paraqus) paths need to be set
at the top of each script.

## Quick start

The following minimal workflow builds a beam unit cell, homogenizes it
under a uniaxial macroscopic strain and reads off a column of the
effective stiffness. It has to be run inside Abaqus
(`abaqus cae noGUI=quickstart.py`).

```python
import numpy as np
import gmsh

from abqhom.gmshmodel import GmshRVE
from abqhom.abqmodel import AbaqusRVE, LoadHistory
from abqhom.utils import voigt_to_tensor, tensor_to_voigt

# 1. Build the unit cell in Gmsh with a physical group for the struts
gmsh.initialize()
gmsh.model.add("CELL")
# ... add points and lines here ...
group = gmsh.model.add_physical_group(1, [...])
gmsh.model.set_physical_name(1, group, "STRUTS")

# 2. Prepare the RVE (boundary and periodic mesh detected automatically)
rve = GmshRVE("CELL")
rve.add_periodic_mesh_constraints()
rve.generate_mesh(element_size=0.1)

# 3. Import into Abaqus and assign a material and a beam section
unit_cell = AbaqusRVE(rve)
material = unit_cell.model.Material(name="STEEL")
material.Elastic(table=[[210.0, 0.3]])
unit_cell.assign_beam_section("STRUTS", "STEEL", radius=0.05)

# 4. Prescribe a macroscopic strain (here eps_xx) and run the analysis
def strain(time):
    eps = np.zeros(6)
    eps[0] = 1.0e-3
    return time*voigt_to_tensor(eps, "STRAIN")

history = LoadHistory([0.0, 1.0], strain)
unit_cell.model.StaticStep(name="LOAD", previous="Initial")
unit_cell.apply_periodic_BC("LOAD", history)
unit_cell.run_analysis(field_outputs=["U"])

# 5. Effective stress -> a column of the effective stiffness matrix
stress = tensor_to_voigt(unit_cell.average_stress("LOAD"), "STRESS",
                         check_symmetry=False)
print(stress/1.0e-3)
```

## Examples

Ready-to-run examples live in [`examples/`](examples) (see
[`examples/README.md`](examples/README.md)); run each from its own folder
with `abaqus cae noGUI=main.py`:

| Example | RVE | Check |
|---|---|---|
| `lattice_simple_cubic` | simple-cubic beam lattice | analytical `C11 = E·π·r²/L²` |
| `lattice_bcc` | cubic + body-diagonal struts | cubic symmetry |
| `lattice_honeycomb_2d` | regular honeycomb | in-plane isotropy |
| `lattice_kelvin_3d` | Kelvin (truncated octahedron) | cubic symmetry |
| `lattice_plasticity` | elasto-plastic beam lattice | nonlinear response |
| `fiber_reinforced_2d` | fibers in a matrix | Voigt/Reuss bounds |
| `spherical_inclusion_3d` | spherical inclusion (manual boundary) | Voigt/Reuss bounds |
| `beam_fibers_in_matrix_3d` | beam fibers embedded in a continuum | transverse isotropy |
| `gyroid_shell_microgen_3d` | sheet-gyroid TPMS (shell), meshed with microgen | cubic symmetry |

## Tests

Unit tests are in `src/abqhom/tests`, split into a Gmsh-only group
(`tests_common`) and a group that needs the Abaqus kernel
(`tests_abaqus`).

- Everything, from the Abaqus Python interpreter:

  ```python
  import abqhom.tests as tests
  tests.run_abaqus_tests()
  ```

- The Gmsh-only subset, from any Python interpreter with Gmsh:

  ```python
  import abqhom.tests as tests
  tests.run_python_tests()
  ```

Both write a report to the current working directory. If Gmsh is not on
the path, set the environment variable `GMSH_LIB` to the folder
containing `gmsh.py`.

## Repository layout

```
src/abqhom/
    gmshmodel.py    GmshRVE, Boundary  -- geometry, boundary, periodic mesh
    abqmodel.py     AbaqusRVE, LoadHistory  -- model, BCs, analysis, averaging
    utils.py        Voigt/tensor conversion, Voigt/Reuss bounds, I/O helpers
    constants.py    boundary-group definitions and element-type maps
    tests/          unit tests (tests_common, tests_abaqus)
examples/           self-contained example scripts
```

## Author

Jonathan Stollberg

## License

Distributed under the terms of the GNU General Public License v3.0. See
[`LICENSE`](LICENSE) for details.

## Acknowledgements

`abqhom` builds on [Gmsh](https://gmsh.info/) for meshing and on
[Paraqus](https://github.com/InstituteOfMechanics/Paraqus) for the VTK
export.
