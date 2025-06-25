import os
import sys
import random
import numpy as np

base_path = "C:/Users/jonat/Documents/"
gmsh_path = os.path.join(base_path, "gmsh", "lib")
abqhom_path = os.path.join(base_path, "abqhom", "src")
paraqus_path = os.path.join(base_path, "Paraqus", "src")
for i in (gmsh_path, abqhom_path, paraqus_path):
    sys.path.append(i)

# Import Gmsh to create the geometry
import gmsh

# Import homogenization tools
from abqhom.gmshmodel import GmshRVE, Boundary
from abqhom.abqmodel import AbaqusRVE, LoadHistory
from abqhom.utils import reuss, voigt, export_to_csv_file


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
# Set working directory
workdir = os.path.join(base_path,
                       "abqhom/examples/spherical_inclusion_3d/abq")
if not os.path.isdir(workdir):
    os.mkdir(workdir)
os.chdir(workdir)

# Create folders to store results
result_path = os.path.join(workdir, "..", "results")
if not os.path.isdir(result_path):
    os.makedirs(os.path.join(result_path, "reuss"))
    os.makedirs(os.path.join(result_path, "voigt"))
    os.makedirs(os.path.join(result_path, "numerical"))

debug_geometry = False

# RVE size
dx = dy = dz = 10.0
inclusion_fraction = 0.3  # Change the volume fraction of the inclusion here

# Boundary conditions
boundary_conditions = "PERIODIC"  # PERIODIC or DIRICHLET

# Material parameters: Young's modulus - Poisson's ratio
parameters = {"MATRIX":     [3.5, 0.35],   # Epoxy matrix
              "INCLUSION":  [230.0, 0.2]}  # Carbon inclusion

# Name folder to store VTK files
vtk_path = os.path.join(workdir, "..", "vtk_output", str(inclusion_fraction))

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
        return time*strain

    history = LoadHistory(np.array([0.0, 1.0]), strain_history)
    print("Strain = " + str(history(1.0).round(3).tolist()))

    # Create geometry in Gmsh and dump it to an input file
    model_name = "RVE_{}".format(i + 1)
    boundary = create_geometry(model_name, dx=dx, dy=dy, dz=dz,
                               inclusion_fraction=inclusion_fraction)
    unit_cell = GmshRVE(model_name, Boundary(model_name, boundary,
                                             enforce_periodic=True))
    unit_cell.generate_mesh(element_size=1.25)
    input_file = unit_cell.write_abaqus_input()

    if debug_geometry:
        gmsh.fltk.run()
        break

    # Load the model into Abaqus
    unit_cell = AbaqusRVE(unit_cell, input_file)

    # Assign material
    for name in parameters:
        E, nu = parameters[name]
        material = unit_cell.model.Material(name=name)
        material.Elastic(table=[[E, nu]])
        unit_cell.assign_solid_section(name, name, "ELSET_" + name)

    # Define load step
    step_name = "LOAD"
    unit_cell.model.StaticStep(name=step_name, previous="Initial")

    # Apply boundary conditions and run analysis
    if boundary_conditions.upper() == "PERIODIC":
        unit_cell.apply_periodic_BC(step_name, history)
    elif boundary_conditions.upper() == "DIRICHLET":
        unit_cell.apply_affine_displacement_BC(step_name, history)
    else:
        raise RuntimeError("Invalid boundary conditions.")

    unit_cell.run_analysis(field_outputs=["NFORC", "U", "S"])
    stress = unit_cell.compute_averaged_stress(step_name)
    print("Stress = " + str(stress.round(3).tolist()))

    # Compute effective stiffness tensor
    C[:, i] = stress/strain_magnitude

    # Export to VTK
    unit_cell.export_to_vtk(vtk_path, field_outputs=["U", ("S", "MISES")],
                            sets=["ELSET_MATRIX", "ELSET_INCLUSION"])

gmsh.finalize()

print("\n------------------------------")
print("Effective constitutive matrix:")
print(C.round(3))

# Compute Reuss and Voigt bounds
C_matrix = isotropic_stiffness_tensor(*parameters["MATRIX"])
C_inclusion = isotropic_stiffness_tensor(*parameters["INCLUSION"])
C_reuss = reuss(inclusion_fraction, C_inclusion, C_matrix)
C_voigt = voigt(inclusion_fraction, C_inclusion, C_matrix)

# Compute eigenvalues for comparison
eigen_reuss = np.sort(np.linalg.eig(C_reuss)[0])
eigen_voigt = np.sort(np.linalg.eig(C_voigt)[0])
eigen_numerical = np.sort(np.linalg.eig(C)[0])

# Export eigenvalues
file = str(inclusion_fraction) + ".csv"
export_to_csv_file(eigen_reuss, os.path.join(result_path, "reuss", file))
export_to_csv_file(eigen_voigt, os.path.join(result_path, "voigt", file))
export_to_csv_file(eigen_numerical, os.path.join(result_path, "numerical",
                                                 file))

print("\n------------------------------")
print("Eigenvalues:")
print("Voigt upper bound:     " + str(eigen_voigt.round(3).tolist()))
print("Numerical:             " + str(eigen_numerical.round(3).tolist()))
print("Reuss lower bound:     " + str(eigen_reuss.round(3).tolist()))
