# Manuscript-reproduction examples

These examples reproduce the numerical studies in the accompanying
computational-homogenization manuscript. In contrast to the concise examples
in the parent directory, they contain analysis and verification scripts for
complete multi-stage workflows.

| Directory | Study | Main capabilities |
|---|---|---|
| [`01_small_strain_verification`](01_small_strain_verification) | Small-strain verification | Mixed continuum–shell–beam cell, ownership reduction, analytical stiffness, and Hill–Mandel check |
| [`02_hyperelastic_octet_lattice`](02_hyperelastic_octet_lattice) | Hyperelastic Octet lattice | Periodic buckling mode, eigenmode imperfection, finite-strain response, contact-clearance check, and supercell comparison |
| [`03_elastoplastic_sandwich_cell`](03_elastoplastic_sandwich_cell) | Elastoplastic sandwich cell | Cyclic response and initial/evolved effective yield surfaces |

The scripts use the same `PYTHONPATH` setup for Gmsh and optional Paraqus as
the tutorials. Each example directory contains a README that describes its
execution order, generated results, and any additional dependencies.
