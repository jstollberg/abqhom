# abqhom

[![License: GPL v3+](https://img.shields.io/badge/License-GPL_v3%2B-blue.svg)](LICENSE)
[![Documentation Status](https://readthedocs.org/projects/abqhom/badge/?version=latest)](https://abqhom.readthedocs.io/en/latest/?badge=latest)

`abqhom` is an open-source Python framework for first-order computational
homogenization of periodic two- and three-dimensional unit cells. It connects
geometry and periodic meshing in [Gmsh](https://gmsh.info/) with finite element
analysis in [Abaqus](https://www.3ds.com/products/simulia/abaqus). Continuum,
truss, beam, and shell constituents may occur separately or together in
mixed-dimensional cells.

The package provides the scale-transition machinery around an Abaqus
analysis. Abaqus remains responsible for element formulations, constitutive
integration, nonlinear solution procedures, and output databases.

## Capabilities

- Periodic and affine boundary conditions for two- and three-dimensional
  cells.
- Linear, finite-strain, nonlinear-elastic, and path-dependent analyses using
  Abaqus material models and solution procedures.
- Periodic translational and rotational constraints for continuum and
  structural finite elements.
- Automatic removal of redundant structural-element copies on secondary
  periodic boundaries without constitutive-parameter scaling.
- Effective first Piola–Kirchhoff stress recovery from element nodal forces
  and incremental microscopic work for discrete Hill–Mandel checks.
- Nonconforming structural–continuum coupling through Abaqus embedded
  constraints.
- VTK export through
  [Paraqus](https://github.com/InstituteOfMechanics/Paraqus).

## Scope

Two-dimensional cells must be rectangular in the global X1–X2 plane with
X3 = 0. Three-dimensional cells must be cuboidal and aligned with the global
coordinate axes. Periodic meshes require one-to-one matching nodes on
opposite boundaries. The effective medium is a first-order Cauchy continuum,
and imposed macroscopic curvature and element warping degrees of freedom are
not supported.

See the complete
[modeling assumptions and limitations](https://abqhom.readthedocs.io/en/latest/limitations.html)
before applying the framework to a new cell topology.

## Requirements

- Python 3.9 or newer for package installation and runtime creation.
- Abaqus 2018 or newer with Abaqus/CAE and the solver.
- Gmsh 4.9.4 or newer.
- NumPy as supplied with Abaqus.
- Paraqus 1.0.2 or newer.
- `backports.weakref` for legacy Abaqus releases using Python 2.7.

The source supports both legacy Python 2.7-based Abaqus releases and current
Python 3-based releases. Larger examples may require a full Abaqus license.

## Installation

Until the first PyPI release, install from a source checkout and create the
isolated Abaqus runtime with the supplied helper:

```console
git clone https://github.com/jstollberg/abqhom.git
cd abqhom
python -m pip install .
abqhom-install-abaqus --abqhom-source .
```

The helper does not modify the Abaqus installation. It creates a versioned
runtime below `~/.abqhom/abaqus-runtime` containing only `abqhom`, Gmsh,
Paraqus, and the weak-reference backport. Activate the generated script in the
shell from which Abaqus will be launched.

The [installation guide](https://abqhom.readthedocs.io/en/latest/installation.html)
describes activation on Windows and Linux, cluster usage, and setup
verification.

## Quick start

After activating the runtime, run the simple-cubic beam-lattice tutorial from
its directory:

```console
cd examples/beam_lattice_simple_cubic_3d
abaqus cae noGUI=main.py
```

The script constructs and periodically meshes the lattice, evaluates six
independent small-strain cases, compares the axial stiffness with its
analytical value, and writes CSV and VTK output.

See the
[`beam_lattice_simple_cubic_3d/main.py`](examples/beam_lattice_simple_cubic_3d/main.py)
source for the minimal complete workflow and the
[`examples` overview](examples/README.md) for the full tutorial collection.

## Examples

The concise tutorials under [`examples/`](examples) cover beam and shell
lattices, continuum composites, embedded mixed-dimensional cells,
finite-strain hyperelasticity, cyclic plasticity, and imported periodic
meshes.

The larger [`paper_examples`](examples/paper_examples) reproduce the numerical
studies from the accompanying publication. Its reference will be added upon
publication.

| Example | Study |
|---|---|
| [`01_small_strain_verification`](examples/paper_examples/01_small_strain_verification) | Mixed continuum–shell–beam verification |
| [`02_hyperelastic_octet_lattice`](examples/paper_examples/02_hyperelastic_octet_lattice) | Imperfect hyperelastic Octet lattice and supercell comparison |
| [`03_elastoplastic_sandwich_cell`](examples/paper_examples/03_elastoplastic_sandwich_cell) | Cyclic plasticity and effective yield surfaces |

Execution order and generated data products are documented in the
[`paper_examples` README](examples/paper_examples/README.md).

## Documentation and support

The [online documentation](https://abqhom.readthedocs.io/en/latest/) contains
the installation guide, modeling assumptions, user workflow, examples, API
reference, and testing instructions. Questions, bug reports, and feature
requests can be submitted through
[GitHub Issues](https://github.com/jstollberg/abqhom/issues).
Contributions are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md) for the
development workflow and submission guidelines.

## Citation

If you use `abqhom` in your research, please cite the following publication:

```bibtex
% Reference will be added upon publication.
```

## Development disclosure

Generative AI tools assisted with code review, refactoring, documentation,
and language editing. All resulting changes were reviewed and tested by the
author, who remains responsible for the software and documentation.

## License

Copyright (C) 2023-2026 Jonathan Stollberg.

Distributed under the terms of the GNU General Public License v3.0 or later.
See [`LICENSE`](LICENSE) for the license text.
