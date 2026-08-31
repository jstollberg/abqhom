"""
Homogenization of a 3d two-phase composite with a spherical inclusion.

The effective stiffness is bracketed by the Voigt (upper) and Reuss
(lower) bounds, computed here for comparison. Because the automatic
boundary detection fails for this OpenCASCADE geometry (the sphere cut
breaks the point-in-entity test), the RVE boundary is provided manually
through the optional 'boundary' argument of GmshRVE -- an example of how
to fall back to a user-defined boundary.

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

from abqhom.utils import (tensor_to_voigt, voigt_to_tensor, reuss, voigt,
                          export_to_csv_file)
from abqhom.gmshmodel import GmshRVE, Boundary
from abqhom.abqmodel import AbaqusRVE, LoadHistory


def create_geometry(model_name, dx=1.0, dy=1.0, dz=1.0,
                    inclusion_fraction=0.2):
    # Compute the radius of the inclusion
    radius = ((0.75*inclusion_fraction*dx*dy*dz)/np.pi)**(1.0/3.0)
    if 2*radius >= min(dx, dy, dz):
        raise ValueError("Inclusion is bigger than the RVE.")

    # Initialize Gmsh and add a new model
    if not gmsh.is_initialized():
        gmsh.initialize()
    gmsh.model.add(model_name)

    gmsh.option.set_number("Geometry.Tolerance", 1e-5)

    # Center point should be coordinate system origin
    x = y = z = 0.0

    # Create the model
    gmsh.model.occ.add_box(x - dx/2.0, y - dy/2.0, z - dz/2.0, dx, dy, dz)
    gmsh.model.occ.add_sphere(x, y, z, radius)
    gmsh.model.occ.cut([(3, 1)], [(3, 2)])
    gmsh.model.occ.add_surface_loop([7])
    gmsh.model.occ.add_volume([1])
    gmsh.model.occ.synchronize()

    # Create map containing information on mesh periodicity and physical groups
    # This is needed since automated generation is not possible for this
    # example
    boundary = {"POINT_1": [(0, 12)],
                "POINT_2": [(0, 15)],
                "POINT_3": [(0, 18)],
                "POINT_4": [(0, 14)],
                "POINT_5": [(0, 11)],
                "POINT_6": [(0, 16)],
                "POINT_7": [(0, 17)],
                "POINT_8": [(0, 13)],
                "EDGE_1":  [(1, 20)],
                "EDGE_2":  [(1, 27)],
                "EDGE_3":  [(1, 25)],
                "EDGE_4":  [(1, 19)],
                "EDGE_5":  [(1, 22)],
                "EDGE_6":  [(1, 24)],
                "EDGE_7":  [(1, 23)],
                "EDGE_8":  [(1, 17)],
                "EDGE_9":  [(1, 16)],
                "EDGE_10": [(1, 21)],
                "EDGE_11": [(1, 26)],
                "EDGE_12": [(1, 18)],
                "FACE_1":  [(2, 12)],
                "FACE_2":  [(2, 10)],
                "FACE_3":  [(2, 9)],
                "FACE_4":  [(2, 11)],
                "FACE_5":  [(2, 8)],
                "FACE_6":  [(2, 13)]
                }

    # Add physical groups
    matrix = gmsh.model.add_physical_group(3, [1])
    gmsh.model.set_physical_name(3, matrix, "MATRIX")
    matrix = gmsh.model.add_physical_group(3, [2])
    gmsh.model.set_physical_name(3, matrix, "INCLUSION")

    return boundary


def isotropic_stiffness_tensor(E, nu):
    C = np.array([[1 - nu, nu, nu, 0.0, 0.0, 0.0],
                  [nu, 1 - nu, nu, 0.0, 0.0, 0.0],
                  [nu, nu, 1 - nu, 0.0, 0.0, 0.0],
                  [0.0, 0.0, 0.0, 0.5*(1 - 2*nu), 0.0, 0.0],
                  [0.0, 0.0, 0.0, 0.0, 0.5*(1 - 2*nu), 0.0],
                  [0.0, 0.0, 0.0, 0.0, 0.0, 0.5*(1 - 2*nu)]])
    C *= E/((1 + nu)*(1 - 2*nu))

    return C


# ---------------------------------------------------------------
# Working directory (created next to this script)
workdir = os.path.join(here, "abq")
if not os.path.isdir(workdir):
    os.makedirs(workdir)
os.chdir(workdir)

# RVE and material parameters (E - nu)
dx = dy = dz = 10.0
inclusion_fraction = 0.3       # volume fraction of the inclusion
parameters = {"MATRIX":    [3.5, 0.35],    # epoxy matrix
              "INCLUSION": [230.0, 0.2]}   # carbon inclusion

# ---------------------------------------------------------------
# Run homogenization routine
C = np.empty((6, 6))
for i in range(6):
    print("\n------------------------------")
    print("Run case {}.".format(i + 1))

    # Define load history
    strain_magnitude = 0.1 if i < 3 else 0.2

    def strain_history(time):
        strain = np.zeros(6)
        strain[i] = strain_magnitude
        strain = voigt_to_tensor(strain, "STRAIN")
        return time*strain

    history = LoadHistory(np.array([0.0, 1.0]), strain_history)
    print("Strain = " + str(tensor_to_voigt(history(1.0), "STRAIN")
                            .round(3).tolist()))

    # Create geometry in Gmsh and dump it to an input file
    model_name = "RVE_{}".format(i + 1)
    boundary = create_geometry(model_name, dx=dx, dy=dy, dz=dz,
                               inclusion_fraction=inclusion_fraction)
    unit_cell = GmshRVE(model_name, Boundary(model_name, boundary))
    unit_cell.add_periodic_mesh_constraints()
    unit_cell.generate_mesh(element_size=2.0)

    # Load the model into Abaqus and assign the phase materials
    unit_cell = AbaqusRVE(unit_cell)
    for name, (E, nu) in parameters.items():
        material = unit_cell.model.Material(name=name)
        material.Elastic(table=[[E, nu]])
        unit_cell.assign_solid_section(name, name)

    unit_cell.model.StaticStep(name="LOAD", previous="Initial")
    unit_cell.apply_periodic_BC("LOAD", history)
    unit_cell.run_analysis(field_outputs=["U", "S"])

    stress = tensor_to_voigt(unit_cell.average_stress("LOAD"), "STRESS",
                             check_symmetry=False)
    C[:, i] = stress/strain_magnitude

    unit_cell.export_to_vtk(os.path.join(here, "vtk_output"),
                            field_outputs=["U", ("S", "MISES")],
                            element_sets=["MATRIX", "INCLUSION"], clear=(i == 0))
    gmsh.finalize()

# Voigt / Reuss bounds of the two-phase material
C_matrix = isotropic_stiffness_tensor(*parameters["MATRIX"])
C_inclusion = isotropic_stiffness_tensor(*parameters["INCLUSION"])
C_reuss = reuss(inclusion_fraction, C_inclusion, C_matrix)
C_voigt = voigt(inclusion_fraction, C_inclusion, C_matrix)

eig = lambda M: np.sort(np.linalg.eigvals(M)).real.round(3).tolist()
print("\n------------------------------")
print("Effective constitutive matrix (Voigt):")
print(np.array2string(C, precision=4, suppress_small=True))
print("\nEigenvalues (numerical must lie between the bounds):")
print("  Voigt upper bound : {}".format(eig(C_voigt)))
print("  numerical         : {}".format(eig(C)))
print("  Reuss lower bound : {}".format(eig(C_reuss)))

# Store the effective stiffness
result_dir = os.path.join(here, "results")
if not os.path.isdir(result_dir):
    os.makedirs(result_dir)
export_to_csv_file(C, os.path.join(result_dir, "stiffness.csv"))
