# abqhom
Simple toolbox for linear computational homogenization using Abaqus.

Gmsh and the SDK are needed: https://gmsh.info/. See examples for information on how to use this package.

## Tests

Unit tests live in `src/abqhom/tests` and are split into a Gmsh-only group
(`tests_common`) and a group that needs the Abaqus kernel (`tests_abaqus`).

- Everything, from the Abaqus Python interpreter:
  `import abqhom.tests as tests; tests.run_abaqus_tests()`
- The Gmsh-only subset, from any Python with Gmsh available:
  `import abqhom.tests as tests; tests.run_python_tests()`

Both write a report to the current working directory. If Gmsh is not on the
path, set the environment variable `GMSH_LIB` to the folder containing
`gmsh.py`.