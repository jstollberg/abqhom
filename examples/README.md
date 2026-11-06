# abqhom examples

The examples are divided into concise user tutorials and the more extensive
manuscript-reproduction cases in `paper_examples`. Each tutorial is a
self-contained `main.py`. Run it from its own folder with

```console
abaqus cae noGUI=main.py
```

Before launching Abaqus, activate the isolated runtime created with
`abqhom-install-abaqus`. The examples import abqhom and its dependencies from
that runtime. Generated files are kept next to the script in `abq/` (Abaqus
job files), `results/` (CSV data and, where applicable, plots) and, when VTK
export is part of the example, `vtk_output/`.

All examples use `RVE` as the internal Gmsh and Abaqus model name. Descriptive
names are reserved for constituents, load cases and result directories.

See the [installation guide](https://abqhom.readthedocs.io/en/latest/installation.html)
for runtime creation and activation on Windows and Linux.

## Beam lattices

| Example | Cell | Demonstrates / exercises |
|---|---|---|
| [`beam_lattice_simple_cubic_3d`](beam_lattice_simple_cubic_3d) | edge struts of the cube | retaining the complete periodic cell through linear-elastic ownership weighting, periodic beam constraints and analytical axial-stiffness verification |
| [`beam_lattice_honeycomb_2d`](beam_lattice_honeycomb_2d) | regular hexagonal honeycomb | 2D boundary classification and pairing for struts intersecting cell edges away from the corners, and in-plane isotropy |
| [`beam_lattice_kelvin_3d`](beam_lattice_kelvin_3d) | Kelvin (truncated octahedron) | finite-strain beam homogenization for a cell whose vertices lie on faces rather than only at corners |
| [`beam_lattice_simple_cubic_plasticity_3d`](beam_lattice_simple_cubic_plasticity_3d) | simple cubic, cyclic perfect plasticity | path-dependent beam response, cyclic loading and microscopic plastic-strain export |

The simple-cubic example retains all periodic boundary members and scales
their Young's modulus by the edge ownership factor (1/4). For linear
elasticity, this weights the complete beam response. The Kelvin cell instead
uses the automatic removal of redundant secondary element copies.

## Continuum and mixed-dimensional cells

| Example | Cell | Demonstrates / exercises |
|---|---|---|
| [`continuum_fiber_composite_2d`](continuum_fiber_composite_2d) | fibers in a matrix (plane strain) | affine displacement boundary conditions and verification of the effective stiffness against Voigt and Reuss bounds |
| [`continuum_spherical_inclusion_3d`](continuum_spherical_inclusion_3d) | spherical inclusion | 3D periodic boundary conditions, automatic boundary classification and finite-strain continuum response |
| [`mixed_embedded_beam_fibers_3d`](mixed_embedded_beam_fibers_3d) | beam fibers embedded in a continuum | mixed-dimensional embedded coupling and periodic constraints on embedded boundary nodes |

## Shell cells

| Example | Cell | Demonstrates / exercises |
|---|---|---|
| [`shell_honeycomb_auxetic_3d`](shell_honeycomb_auxetic_3d) | re-entrant honeycomb, shell walls | periodic shell homogenization and the negative in-plane Poisson response of a re-entrant architecture |
| [`shell_gyroid_3d`](shell_gyroid_3d) | sheet-gyroid TPMS with an imported shell mesh | importing a pre-generated periodic mesh and computing its finite-strain shell response |

The gyroid mesh is generated once with the optional
[microgen](https://github.com/3MAH/microgen) package by `generate_mesh.py`.
The Abaqus analysis reads the resulting `gyroid.msh` and does not require
microgen.

## Manuscript-reproduction cases

The three cases under `paper_examples` reproduce the numerical examples of
the accompanying manuscript. They include multi-stage analyses and
verification scripts and are therefore intentionally more detailed than the
tutorials. See the [`paper_examples` README](paper_examples/README.md) for
their execution order and generated output.
