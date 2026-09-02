"""
Homogenization of a simple-cubic beam-lattice unit cell.

The unit cell consists of the twelve struts along the edges of the cubic
bounding box. Because the edge struts are shared by four neighboring
cells, this is the simplest example that exercises the automatic periodic
sharing weights in average_stress.

Analytical check
----------------
Under uniaxial strain only the four struts parallel to the loading
direction carry an axial force N = E*A*eps. Each edge is shared by four
cells (weight 1/4), so their effective contribution is that of a single
strut and

    C11 = E*A / L**2 = E*pi*r**2 / L**2 ,

while the transverse coupling C12 vanishes (no Poisson effect for the
axial struts).

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


def create_geometry(model_name, length):
    """Create the twelve edge struts of a cube of the given edge length."""
    if not gmsh.is_initialized():
        gmsh.initialize()
    if model_name in gmsh.model.list():
        gmsh.model.set_current(model_name)
        gmsh.model.remove()
    gmsh.model.add(model_name)
    g = gmsh.model.geo

    h = 0.5*length
    corners = [g.add_point(x, y, z)
               for x in (-h, h) for y in (-h, h) for z in (-h, h)]

    # Connect the corners that differ in exactly one coordinate (edges)
    struts = []
    coords = [(x, y, z) for x in (-h, h) for y in (-h, h) for z in (-h, h)]
    for i in range(8):
        for j in range(i + 1, 8):
            shared = 0
            for a, b in zip(coords[i], coords[j]):
                if a == b:
                    shared += 1
            if shared == 2:
                struts.append(g.add_line(corners[i], corners[j]))

    g.synchronize()
    group = gmsh.model.add_physical_group(1, struts)
    gmsh.model.set_physical_name(1, group, "STRUTS")


# ---------------------------------------------------------------
# Working directory (created next to this script)
workdir = os.path.join(here, "abq")
if not os.path.isdir(workdir):
    os.makedirs(workdir)
os.chdir(workdir)

# Material and geometry
E = 210.0        # Young's modulus of the strut material
nu = 0.3         # Poisson's ratio
length = 10.0    # edge length of the cubic unit cell
radius = 0.5     # strut radius
element_size = 1.0

# ---------------------------------------------------------------
# Build the model once (geometry, mesh, material, section, step). Only the
# macroscopic strain changes between load cases, so the model is reused and
# only the periodic boundary conditions are updated in the loop.
model_name = "RVE"
create_geometry(model_name, length)
unit_cell = GmshRVE(model_name)
unit_cell.add_periodic_mesh_constraints()
unit_cell.generate_mesh(element_size=element_size)

unit_cell = AbaqusRVE(unit_cell)
material = unit_cell.model.Material(name="STRUT")
material.Elastic(table=[[E, nu]])
unit_cell.assign_beam_section("STRUTS", "STRUT", radius)
unit_cell.model.StaticStep(name="LOAD", previous="Initial")

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
    unit_cell.apply_periodic_BC("LOAD", history)
    unit_cell.run_analysis(field_outputs=["U"])

    stress = tensor_to_voigt(unit_cell.average_stress("LOAD"), "STRESS")
    C[:, i] = stress/strain_magnitude

    unit_cell.export_to_vtk(os.path.join(here, "vtk_output",
                                         "case_{}".format(i + 1)),
                            field_outputs=["U"], clear=True)

gmsh.finalize()

# ---------------------------------------------------------------
C11_analytical = E*np.pi*radius**2/length**2
print("\n------------------------------")
print("Effective constitutive matrix (Voigt):")
print(np.array2string(C, precision=4, suppress_small=True))
print("\nComparison with the analytical simple-cubic result:")
print("  C11 (FE)         = {:.6e}".format(C[0, 0]))
print("  C11 (analytical) = {:.6e}".format(C11_analytical))
print("  relative error   = {:.2e}".format(
    abs(C[0, 0] - C11_analytical)/C11_analytical))
print("  max |C12,C13,C23| (should vanish) = {:.2e}".format(
    max(abs(C[0, 1]), abs(C[0, 2]), abs(C[1, 2]))))

# Store the effective stiffness
result_dir = os.path.join(here, "results")
if not os.path.isdir(result_dir):
    os.makedirs(result_dir)
export_to_csv_file(C, os.path.join(result_dir, "stiffness.csv"))
