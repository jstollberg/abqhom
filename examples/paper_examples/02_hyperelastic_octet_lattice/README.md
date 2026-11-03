# Example 2: hyperelastic Octet lattice with imperfections

This example reproduces the finite-strain Octet-lattice study reported in
Section 5.2 of the accompanying manuscript. The periodic beam lattice has
`r/ell = 0.04`, a compressible Neo-Hookean material, and an eigenmode
imperfection of amplitude `a/r = 0.01`. Compression, tension, and simple
shear are imposed up to a macroscopic deformation of 40%. Beam contact is
not included in the analyses and is checked independently from the stored
deformation histories.

## Execution

Activate the isolated Abaqus runtime as described in the
[installation guide](https://abqhom.readthedocs.io/en/latest/installation.html),
then run the analyses from this directory in the following order:

```console
abaqus cae noGUI=buckling.py
abaqus cae noGUI=main.py
```

The first command computes two periodic buckling modes for the primitive cell
and, by default, for the `2 x 2 x 2` supercell. The first mode defines the
geometric imperfection; the second records the neighboring buckling mode. The
second command applies the translational part of the first mode and runs
compression, tension, and shear for both cell sizes. Set
`include_supercell=False` in the final call of either script to omit the
supercell analysis.

Paraqus generates the VTK histories used by the contact-clearance check.
After exporting the primitive-cell histories, run the check in a Python
environment containing the optional `vtk` Python package. It is not part of
the abqhom Abaqus runtime:

```console
python verify_contact_clearance.py
```

The check compares every retained primitive-cell beam element with the
central cell and its 26 periodic neighbors. Members sharing an idealized
Octet junction are excluded. The script evaluates the exact minimum distance
between all remaining discretized centerlines. It raises an error if contact
occurs or a positive surface clearance cannot be resolved within the adopted
numerical tolerance.

## Supercell check

The shared functions in `octet.py` construct a conforming supercell in which
points and struts on internal cell interfaces occur only once. Periodic
equations and ownership reduction therefore act only on its outer boundary.
Both cell sizes use 16 B31 elements per strut. Their compression, tension,
and shear responses can be compared directly through the respective
`stress_curve.csv` files.

## Generated output

| Directory | Contents |
|---|---|
| `abq/` | Primitive-cell Abaqus jobs |
| `results/` | Primitive-cell modes, responses, and diagnostics |
| `vtk_output/` | Primitive-cell VTK collections |
| `supercell_abq/` | Supercell Abaqus jobs |
| `supercell_results/` | Supercell modes, responses, and diagnostics |
| `supercell_vtk_output/` | Supercell VTK collections |

Each load-case directory contains:

| File | Contents |
|---|---|
| `full_response.csv` | Pseudo-time, loading measure, nine components of the first Piola–Kirchhoff stress, nine components of the Cauchy stress, stress-asymmetry residual, and Hill–Mandel residual |
| `stress_curve.csv` | Loading measure and the corresponding first Piola–Kirchhoff and Cauchy stress components |
| `verification_metrics.csv` | Peak kinetic-energy ratio and its pseudo-time, final kinetic-energy ratio, ratio of the separate maximum kinetic and internal energies, maximum stress-asymmetry and Hill–Mandel residuals, and maximum nodal rotation over the complete history |

Tensor components in `full_response.csv` are stored in row-major order. The
first Hill–Mandel entry is undefined because no preceding frame exists. The
peak kinetic-to-internal-energy ratio excludes the initial interval
`tau < 0.02`. For the primitive cell, `contact_clearance.csv` additionally
contains pseudo-time, loading magnitude, and minimum centerline distance.

Each `mode_XX.csv` contains the node label, three reference coordinates, and
three normalized translational mode components. Primitive-cell
`eigenvalues.csv` stores the mode number, eigenvalue, and outer-periodicity
residual. The supercell file additionally stores the residual with respect to
primitive-cell periodicity and the number of compared node pairs.
