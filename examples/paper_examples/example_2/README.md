# Example 2: finite-strain Octet lattice

The example considers a periodic Octet beam lattice with `r/ell = 0.04`, a
compressible Neo-Hookean material, and an eigenmode imperfection of amplitude
`a/r = 0.01`. Compression, tension, and simple shear are imposed up to a
macroscopic deformation of 40%. Beam contact is not included.

## Execution

Run the files in this order from an Abaqus command prompt:

```text
abaqus cae noGUI=buckling.py
abaqus cae noGUI=main.py
```

The first command generates the primitive-cell modes in
`final_results/buckling` and, by default, the supercell modes below
`supercell_results`. The second command applies the corresponding
imperfections and runs all three nonlinear load cases for both cell sizes.
Set `include_supercell=False` in the final call of either script to omit the
supercell analysis. Primitive-cell Abaqus files, tabulated results, and VTK
collections are written to `final_abq`, `final_results`, and `final_vtk`.
The corresponding supercell outputs are written to `sc_abq`,
`supercell_results`, and `supercell_vtk`.

Paraqus is optional for the analyses but required to generate the VTK files
used by the clearance check. Run the check in a Python environment containing
VTK after the primitive-cell VTK histories have been exported:

```text
python verify_contact_clearance.py
```

The script checks every retained primitive-cell beam against the central cell
and its 26 periodic neighbors. It excludes members sharing an idealized Octet
junction, determines the exact minimum distance between all remaining
discretized centerlines, and raises an error if the corresponding surface
clearance is not positive.

## Supercell check

The shared utilities in `octet.py` construct a conforming supercell in which
points and struts on internal cell interfaces occur only once. The periodic
equations and ownership reduction therefore act only on the outer supercell
boundary. Both cell sizes use 16 B31 elements per strut.
Buckling modes and the compression, tension, and shear responses are written
below `supercell_results`.

The common discretization permits a direct cell-size comparison through the
respective `stress_curve.csv` files.

## Generated output

| Directory | Contents |
|---|---|
| `final_abq/` | Primitive-cell Abaqus jobs |
| `final_results/` | Primitive-cell modes, responses, and diagnostics |
| `final_vtk/` | Primitive-cell VTK collections |
| `sc_abq/` | Supercell Abaqus jobs |
| `supercell_results/` | Supercell modes, responses, and diagnostics |
| `supercell_vtk/` | Supercell VTK collections |

Each load-case directory contains:

| File | Contents |
|---|---|
| `full_response.csv` | Pseudo-time, loading measure, nine components of `P`, nine components of `sigma`, stress-asymmetry residual, and Hill--Mandel residual |
| `stress_curve.csv` | Loading measure and the corresponding components of `P` and `sigma` |
| `verification_metrics.csv` | Peak kinetic-energy ratio and its pseudo-time, final and global kinetic-energy ratios, maximum stress-asymmetry and Hill--Mandel residuals, and maximum nodal rotation |

Tensor components in `full_response.csv` are stored in row-major order. The
first Hill--Mandel entry is undefined because no preceding frame exists.
For the primitive cell, `contact_clearance.csv` additionally contains
pseudo-time, loading magnitude, and the minimum centerline distance.

Each `mode_XX.csv` contains the node label, three reference coordinates, and
three normalized translational mode components. Primitive-cell
`eigenvalues.csv` stores the mode number, eigenvalue, and periodicity
residual. The supercell file additionally stores the residual with respect to
primitive-cell periodicity and the number of compared node pairs.

All Abaqus files, tabulated results, and VTK collections are generated on
demand and are not retained in the repository.
