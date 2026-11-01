# Example 1: mixed-dimensional small-strain verification

This example verifies the homogenization framework on a conforming periodic
cell containing continuum, shell, and beam elements. It reproduces the study
reported in Section 5.1 of the accompanying manuscript.

## Verification targets

The analysis evaluates:

- the complete effective stiffness from six independent strain cases;
- the analytical normal-stiffness block;
- major and cubic symmetry residuals;
- the discrete Hill–Mandel residual in every strain case;
- equivalence between the reduced cell and the complete cell with
  ownership-weighted linear-elastic properties.

Homogenization uses the ownership-consistent reduced cell. The complete
periodic cell enters only the independent ownership-weighted comparison.

## Execution

Expose Gmsh and optional Paraqus through `PYTHONPATH` as described in the
project documentation. Run the analysis from this directory:

```text
abaqus cae noGUI=main.py
```

Paraqus is required only for VTK export. If it is unavailable, the numerical
verification and CSV export still run.

## Generated output

| Directory | Contents |
|---|---|
| `abq/` | Reduced-cell Abaqus jobs |
| `abq_weighted/` | Complete ownership-weighted comparison jobs |
| `results/` | Stiffness matrices and verification residuals |
| `vtk_output/` | Reduced and ownership-weighted VTK collections for all six strain cases |

The `results/` directory contains:

| File | Contents |
|---|---|
| `effective_stiffness_computed.csv` | Numerical 6 x 6 stiffness matrix |
| `normal_stiffness_computed.csv` | Numerical 3 x 3 normal block |
| `normal_stiffness_analytical.csv` | Analytical 3 x 3 normal block |
| `hill_mandel_residuals.csv` | Hill–Mandel residuals for the six strain cases |
| `effective_stiffness_weighted.csv` | Stiffness of the complete ownership-weighted cell |
| `verification_residuals.csv` | Summary residuals reported for the verification study |

Rows and columns of the stiffness matrices use the Voigt ordering
`(11, 22, 33, 23, 13, 12)`. Entries of `hill_mandel_residuals.csv` follow the
same strain-case order.

The columns of `verification_residuals.csv` contain the normal-block error,
major-symmetry residual, cubic-symmetry residual, reduction residual, and
maximum Hill–Mandel residual. The CSV header records these names explicitly.

All model and material parameters are defined together near the beginning of
`main.py`. The analysis uses 25 elements per cell edge and a macroscopic
strain magnitude of 0.01.
