"""
Linear homogenization of a 3d Kelvin (tetrakaidecahedron) beam-lattice
unit cell.

Like the honeycomb example, this demonstrates a general lattice cell
whose struts do not connect to the corners of the bounding box: the
Kelvin cell's vertices lie on the faces of the cubic bounding box. The
workflow is identical to the classic "cubic" cells -- the boundary, its
periodic connections and the rigid-body anchor are found automatically.

Run with:  abaqus cae noGUI=main.py
"""

import os
import sys
import inspect
from itertools import permutations, product

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

from abqhom.utils import voigt_to_tensor, tensor_to_voigt, export_to_csv_file
from abqhom.gmshmodel import GmshRVE
from abqhom.abqmodel import AbaqusRVE, LoadHistory


def create_kelvin(model_name, edge=4.0):
    """Create a Kelvin (truncated octahedron) beam lattice in the cube
    [-edge/2, edge/2]^3. The 24 vertices are the permutations of
    (0, +-1, +-2)*edge/4; each has a +-2 coordinate, hence lies on a cube
    face. Struts connect vertices at distance sqrt(2)*edge/4.
    """
    if not gmsh.is_initialized():
        gmsh.initialize()
    if model_name in gmsh.model.list():
        gmsh.model.set_current(model_name)
        gmsh.model.remove()
    gmsh.model.add(model_name)
    g = gmsh.model.geo

    s = edge/4.0
    verts = set()
    for perm in set(permutations((0, 1, 2))):
        for signs in product((1, -1), repeat=3):
            verts.add(tuple(round(perm[i]*signs[i]*s, 9) for i in range(3)))
    verts = sorted(verts)

    tag = {v: g.add_point(v[0], v[1], v[2]) for v in verts}
    target = np.sqrt(2.0)*s
    lines = []
    for i in range(len(verts)):
        for j in range(i + 1, len(verts)):
            d = np.linalg.norm(np.array(verts[i]) - np.array(verts[j]))
            if abs(d - target) < 1e-6:
                lines.append(g.add_line(tag[verts[i]], tag[verts[j]]))

    g.synchronize()
    grp = gmsh.model.add_physical_group(1, lines)
    gmsh.model.set_physical_name(1, grp, "STRUTS")


# ---------------------------------------------------------------
# Working directory (created next to this script)
workdir = os.path.join(here, "abq")
if not os.path.isdir(workdir):
    os.makedirs(workdir)
os.chdir(workdir)

# Material parameters (of the solid strut material)
E = 1.0     # Young's modulus
nu = 0.3    # Poisson's ratio

# Lattice parameters
edge = 4.0          # cubic unit-cell edge length
radius = 0.05       # strut cross-section radius
element_size = 0.6  # beam element length

# ---------------------------------------------------------------
# Run the homogenization routine (6 load cases in 3d)
components = [0, 1, 2, 3, 4, 5]  # Voigt: xx, yy, zz, yz, xz, xy
strain_magnitude = 1.0e-3

C = np.empty((6, 6))
for col, comp in enumerate(components):
    print("\n------------------------------")
    print("Run load case {}.".format(col + 1))

    def strain_history(time):
        strain = np.zeros(6)
        strain[comp] = strain_magnitude
        return time*voigt_to_tensor(strain, "STRAIN")

    history = LoadHistory(np.array([0.0, 1.0]), strain_history)

    # Create geometry and periodic mesh in Gmsh
    model_name = "RVE_{}".format(col + 1)
    create_kelvin(model_name, edge=edge)
    unit_cell = GmshRVE(model_name)
    unit_cell.add_periodic_mesh_constraints()
    unit_cell.generate_mesh(element_size=element_size)

    # Load the model into Abaqus and assign a beam section to the struts
    unit_cell = AbaqusRVE(unit_cell)
    material = unit_cell.model.Material(name="LINEAR_ELASTICITY")
    material.Elastic(table=[[E, nu]])
    unit_cell.assign_beam_section("STRUTS", "LINEAR_ELASTICITY", radius)

    # Add load step and periodic boundary conditions (anchor is picked
    # automatically -- no corner point required)
    unit_cell.model.StaticStep(name="LOAD", previous="Initial")
    unit_cell.apply_periodic_BC("LOAD", history)

    # Run the analysis and compute the average stress
    unit_cell.run_analysis(field_outputs=["U"])
    stress = unit_cell.average_stress("LOAD")
    stress = tensor_to_voigt(stress, "STRESS", check_symmetry=False)
    print("Stress = " + str(stress.round(8).tolist()))

    C[:, col] = stress/strain_magnitude

    # Export the deformed struts to VTK (strut radius taken from section)
    unit_cell.export_to_vtk(os.path.join(here, "vtk_output"),
                            field_outputs=["U"], clear=(col == 0))
    gmsh.finalize()

print("\n------------------------------")
print("Effective constitutive matrix (Voigt):")
print(np.array2string(C, precision=6, suppress_small=True))
print("\nChecks for the Kelvin cell (should have cubic symmetry):")
print("  symmetry  max|C - C^T|      = {:.2e}".format(np.max(np.abs(C - C.T))))
print("  C11=C22=C33 spread          = {:.2e}".format(
    np.ptp([C[0, 0], C[1, 1], C[2, 2]])))
print("  C12=C13=C23 spread          = {:.2e}".format(
    np.ptp([C[0, 1], C[0, 2], C[1, 2]])))
print("  C44=C55=C66 spread          = {:.2e}".format(
    np.ptp([C[3, 3], C[4, 4], C[5, 5]])))

# Export the effective stiffness
result_path = os.path.join(here, "results")
if not os.path.isdir(result_path):
    os.makedirs(result_path)
export_to_csv_file(C, os.path.join(result_path, "kelvin_C.csv"))
