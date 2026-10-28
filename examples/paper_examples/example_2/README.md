# Example 2: finite-strain Octet lattice

The example considers a periodic Octet beam lattice with (r/\ell=0.04),
a compressible Neo-Hookean material, and an eigenmode imperfection of
amplitude (a/r=0.01). Compression, tension, and simple shear are imposed up
to a macroscopic deformation of 40%. Beam contact is not included.

Run the files in this order from an Abaqus command prompt:

```text
abaqus cae noGUI=buckling.py
abaqus cae noGUI=main.py
```

The first command generates the periodic eigenmodes in
`final_results/buckling`. The second command applies the imperfection to the
Gmsh mesh and runs all three nonlinear load cases. Abaqus files, tabulated
results, and VTK collections are written to `final_abq`, `final_results`, and
`final_vtk`, respectively.

The response plots can subsequently be created in the `chatgpt` environment:

```text
python plot_results.py
python plot_all_cauchy_components.py
python plot_octet_buckling_modes.py
```

`verify_contact_clearance.py` checks every stored frame against the central
cell and its 26 periodic neighbors. It excludes only members sharing an
idealized Octet junction and verifies the clearance between all remaining
beam centerlines against the diameter (2r).

The repository retains the compact tabulated results and VTK collections
needed for plotting and verification. The much larger Abaqus job directories
are generated on demand and are not retained.

## Supercell check

`supercell.py` constructs a conforming Octet supercell in which points and
struts on internal cell interfaces occur only once. The periodic equations
and ownership reduction therefore act only on the outer supercell boundary.
The production supercell study is run with

```text
abaqus cae noGUI=supercell_buckling.py
abaqus cae noGUI=supercell_compression.py
```

Both scripts use the production resolution of 16 B31 elements per strut by
default. A different common mesh density can be selected through the
`OCTET_ELEMENTS_PER_STRUT` environment variable. Buckling modes and nonlinear
compression curves are written to mesh-tagged directories below
`supercell_results`.

For diagnostic reruns, `OCTET_SUPERCELL_MODE` selects mode 1 or 3 and
`OCTET_STEP_DURATION` controls the duration of the implicit dynamic step.

For a direct cell-size comparison, the single-cell and supercell scripts must
be run with the same `OCTET_ELEMENTS_PER_STRUT` value.
`plot_supercell_comparison.py` compares the production results at 16 elements
per strut.
