"""
Homogenization of a sheet-gyroid shell unit cell, meshed with microgen.

This is a variant of ``gyroid_shell_3d`` that uses the periodic gyroid
surface mesh produced by the microgen package instead of the built-in
marching-tetrahedra routine. The mesh is pre-computed and exported to
``gyroid.msh`` by ``generate_mesh.py`` (run separately in an environment
that has microgen installed), so Abaqus itself does not need microgen -- it
only reads the mesh file.

The pre-built mesh is loaded into Gmsh and handed to abqhom through
``GmshRVE.read_from_mesh``, which classifies the boundary nodes onto the
bounding-box groups by coordinate and checks the mesh for periodicity. The
full effective stiffness (6 load cases) is then homogenized with shell
elements and checked for cubic symmetry.

Run with:  abaqus cae noGUI=main.py
"""

import os
import sys
import inspect

import numpy as np

# ---------------------------------------------------------------
# Locate this example and put the shipped abqhom source first on the
# path. Under 'abaqus cae noGUI=main.py' __file__ is undefined, so the
# script location is obtained via inspect instead.
here = os.path.dirname(os.path.abspath(
    inspect.getframeinfo(inspect.currentframe()).filename))
abqhom_path = os.path.abspath(os.path.join(here, os.pardir, os.pardir, "src"))
gmsh_path = "C:/Users/jonat/Documents/gmsh/lib"        # TODO: set gmsh path
paraqus_path = "C:/Users/jonat/Documents/Paraqus/src"  # optional (VTK export)
for path in (abqhom_path, gmsh_path, paraqus_path):
    if path not in sys.path:
        sys.path.insert(0, path)

import gmsh

from abqhom.utils import (tensor_to_voigt, voigt_to_tensor,
                          export_to_csv_file)
from abqhom.gmshmodel import GmshRVE
from abqhom.abqmodel import AbaqusRVE, LoadHistory

# ---------------------------------------------------------------
# Working directory (created next to this script)
workdir = os.path.join(here, "abq")
if not os.path.isdir(workdir):
    os.makedirs(workdir)
os.chdir(workdir)

# Shell thickness and material
thickness = 0.05     # shell thickness
E = 210.0            # Young's modulus of the shell material
nu = 0.3             # Poisson's ratio

# Load the pre-computed microgen mesh (physical group 'GYROID').
mesh_file = os.path.join(here, "gyroid.msh")
if not os.path.isfile(mesh_file):
    raise RuntimeError(
        "The mesh file '{}' does not exist. Generate it first by running "
        "generate_mesh.py in an environment that has microgen installed."
        .format(mesh_file))

gmsh.initialize()
# Merge the mesh into an explicitly named model (Abaqus uppercases the
# part name, so the model name must match that case).
model_name = "GYROID"
gmsh.model.add(model_name)
gmsh.merge(mesh_file)
unit_cell = GmshRVE.read_from_mesh(model_name)
if not unit_cell.periodic_mesh:
    raise RuntimeError("The imported mesh is not periodic.")

# ---------------------------------------------------------------
# Build the Abaqus model once (mesh, material, section, step). Only the
# macroscopic strain changes between load cases, so the model is reused
# and only the periodic boundary conditions are updated in the loop.
cell = AbaqusRVE(unit_cell)
material = cell.model.Material(name="SHELL")
material.Elastic(table=[[E, nu]])
cell.assign_shell_section("GYROID", "SHELL", thickness)
cell.model.StaticStep(name="LOAD", previous="Initial")

# Run the homogenization routine (6 load cases in 3d)
strain_magnitude = 1.0e-3
C = np.zeros((6, 6))
for i in range(6):
    print("\n------------------------------")
    print("Run load case {}.".format(i + 1))

    def strain_history(time, component=i):
        strain = np.zeros(6)
        strain[component] = strain_magnitude
        return time*voigt_to_tensor(strain, "STRAIN")

    history = LoadHistory(np.array([0.0, 1.0]), strain_history)

    # Re-apply the periodic BC with the new strain (the previous one is
    # cleared automatically) and re-run the analysis on the same model.
    cell.apply_periodic_BC("LOAD", history)
    cell.run_analysis(field_outputs=["U"])

    stress = tensor_to_voigt(cell.average_stress("LOAD"), "STRESS")
    C[:, i] = stress/strain_magnitude

    # All load cases reuse the same RVE (model name 'GYROID'), so the VTK
    # files are written into a separate subfolder per case to avoid
    # overwriting each other.
    cell.export_to_vtk(os.path.join(here, "vtk_output",
                                    "case_{}".format(i + 1)),
                       field_outputs=["U"], clear=True)

gmsh.finalize()

# ---------------------------------------------------------------
# A gyroid is cubic-symmetric: C11 = C22 = C33, C12 = C13 = C23 and
# C44 = C55 = C66. The Zener ratio characterises the anisotropy.
print("\n------------------------------")
print("Effective constitutive matrix (Voigt):")
print(np.array2string(C, precision=4, suppress_small=True))

C11 = np.mean([C[0, 0], C[1, 1], C[2, 2]])
C12 = np.mean([C[0, 1], C[0, 2], C[1, 2]])
C44 = np.mean([C[3, 3], C[4, 4], C[5, 5]])
print("\nCubic-symmetry check (averaged over the equivalent components):")
print("  C11 = C22 = C33 = {:.4f}  (spread {:.2e})".format(
    C11, np.ptp([C[0, 0], C[1, 1], C[2, 2]])))
print("  C12 = C13 = C23 = {:.4f}  (spread {:.2e})".format(
    C12, np.ptp([C[0, 1], C[0, 2], C[1, 2]])))
print("  C44 = C55 = C66 = {:.4f}  (spread {:.2e})".format(
    C44, np.ptp([C[3, 3], C[4, 4], C[5, 5]])))
print("  Zener anisotropy 2*C44/(C11-C12) = {:.4f}".format(
    2.0*C44/(C11 - C12)))
print("  max asymmetry |C-C^T| = {:.2e}".format(np.max(np.abs(C - C.T))))

# Store the effective stiffness
result_dir = os.path.join(here, "results")
if not os.path.isdir(result_dir):
    os.makedirs(result_dir)
export_to_csv_file(C, os.path.join(result_dir, "stiffness.csv"))
