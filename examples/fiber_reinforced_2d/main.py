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
from abqhom.gmshmodel import GmshRVE
from abqhom.abqmodel import AbaqusRVE
from abqhom.utils import reuss, voigt, export_to_csv_file


def create_geometry(model_name, dx=1.0, dy=1.0, fiber_fraction=0.2):
    # Compute fiber radius from volume fraction
    radius = np.sqrt((fiber_fraction*dx*dy)/(2*np.pi))
    if radius >= np.sqrt((dx/2)**2 + (dy/2)**2)/2:
        raise ValueError("Fibers overlap.")

    # Initialize Gmsh and add a new model
    if not gmsh.is_initialized():
        gmsh.initialize()
    gmsh.model.add(model_name)

    gmsh.option.set_number("Geometry.Tolerance", 1e-5)

    # Center point should be coordinate system origin
    x = y = 0.0

    # Create the model
    point_1 = gmsh.model.occ.add_point(x - dx/2, y - dy/2, 0.0)
    point_2 = gmsh.model.occ.add_point(x + dx/2, y - dy/2, 0.0)
    point_3 = gmsh.model.occ.add_point(x + dx/2, y + dy/2, 0.0)
    point_4 = gmsh.model.occ.add_point(x - dx/2, y + dy/2, 0.0)
    point_5 = gmsh.model.occ.add_point(x - dx/2 + radius, y - dy/2, 0.0)
    point_6 = gmsh.model.occ.add_point(x + dx/2 - radius, y - dy/2, 0.0)
    point_7 = gmsh.model.occ.add_point(x + dx/2, y - dy/2 + radius, 0.0)
    point_8 = gmsh.model.occ.add_point(x + dx/2, y + dy/2 - radius, 0.0)
    point_9 = gmsh.model.occ.add_point(x + dx/2 - radius, y + dy/2, 0.0)
    point_10 = gmsh.model.occ.add_point(x - dx/2 + radius, y + dy/2, 0.0)
    point_11 = gmsh.model.occ.add_point(x - dx/2, y + dy/2 - radius, 0.0)
    point_12 = gmsh.model.occ.add_point(x - dx/2, y - dy/2 + radius, 0.0)
    point_13 = gmsh.model.occ.add_point(x - radius, y, 0.0)
    point_14 = gmsh.model.occ.add_point(x + radius, y, 0.0)
    point_15 = gmsh.model.occ.add_point(x, y - radius, 0.0)
    point_16 = gmsh.model.occ.add_point(x, y + radius, 0.0)

    line_1 = gmsh.model.occ.add_line(point_1, point_5)
    line_2 = gmsh.model.occ.add_line(point_5, point_6)
    line_3 = gmsh.model.occ.add_line(point_6, point_2)
    line_4 = gmsh.model.occ.add_line(point_2, point_7)
    line_5 = gmsh.model.occ.add_line(point_7, point_8)
    line_6 = gmsh.model.occ.add_line(point_8, point_3)
    line_7 = gmsh.model.occ.add_line(point_3, point_9)
    line_8 = gmsh.model.occ.add_line(point_9, point_10)
    line_9 = gmsh.model.occ.add_line(point_10, point_4)
    line_10 = gmsh.model.occ.add_line(point_4, point_11)
    line_11 = gmsh.model.occ.add_line(point_11, point_12)
    line_12 = gmsh.model.occ.add_line(point_12, point_1)

    circle_1 = gmsh.model.occ.add_circle_arc(point_5, point_1, point_12)
    circle_2 = gmsh.model.occ.add_circle_arc(point_7, point_2, point_6)
    circle_3 = gmsh.model.occ.add_circle_arc(point_9, point_3, point_8)
    circle_4 = gmsh.model.occ.add_circle_arc(point_11, point_4, point_10)
    circle_5 = gmsh.model.occ.add_circle_arc(point_13, point_15, point_14,
                                             center=False)
    circle_6 = gmsh.model.occ.add_circle_arc(point_14, point_16, point_13,
                                             center=False)

    loop_1 = gmsh.model.occ.add_curve_loop([line_12, line_1, circle_1])
    loop_2 = gmsh.model.occ.add_curve_loop([line_3, line_4, circle_2])
    loop_3 = gmsh.model.occ.add_curve_loop([line_6, line_7, circle_3])
    loop_4 = gmsh.model.occ.add_curve_loop([line_9, line_10, circle_4])
    loop_5 = gmsh.model.occ.add_curve_loop([circle_5, circle_6])
    loop_6 = gmsh.model.occ.add_curve_loop([
        line_2, circle_2, line_5, circle_3,
        line_8, circle_4, line_11, circle_1])

    surface_1 = gmsh.model.occ.add_plane_surface([loop_1])
    surface_2 = gmsh.model.occ.add_plane_surface([loop_2])
    surface_3 = gmsh.model.occ.add_plane_surface([loop_3])
    surface_4 = gmsh.model.occ.add_plane_surface([loop_4])
    surface_5 = gmsh.model.occ.add_plane_surface([loop_5])
    surface_6 = gmsh.model.occ.add_plane_surface([loop_6, loop_5])

    gmsh.model.occ.synchronize()

    # Add physical groups
    matrix = gmsh.model.add_physical_group(2, [surface_6])
    gmsh.model.set_physical_name(2, matrix, "MATRIX")
    fiber = gmsh.model.add_physical_group(2, [surface_1, surface_2, surface_3,
                                              surface_4, surface_5])
    gmsh.model.set_physical_name(2, fiber, "FIBER")

    # Create map containing information on mesh periodicity and physical groups
    boundary_map = {"POINT_1": (0, [point_1]),
                    "POINT_2": (0, [point_2]),
                    "POINT_3": (0, [point_3]),
                    "POINT_4": (0, [point_4]),
                    "EDGE_1":  (1, [line_1, line_2, line_3]),
                    "EDGE_2":  (1, [line_4, line_5, line_6]),
                    "EDGE_3":  (1, [line_9, line_8, line_7]),
                    "EDGE_4":  (1, [line_12, line_11, line_10]),
                    }

    return boundary_map


def isotropic_stiffness_tensor(E, nu):
    C = np.array([[1 - nu, nu, 0.0],
                  [nu, 1 - nu, 0.0],
                  [0.0, 0.0, 0.5*(1 - 2*nu)]])
    C *= E/((1 + nu)*(1 - 2*nu))

    return C


# ---------------------------------------------------------------
# Set working directory
workdir = os.path.join(base_path,
                       "abqhom/examples/fiber_reinforced_2d/abq")
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
dx = dy = 10.0
fiber_fraction = 0.7  # Change the volume fraction of the fibers here

# Boundary conditions
boundary_conditions = "PERIODIC"  # PERIODIC or DIRICHLET

# Material parameters for CFRP: Young's modulus - Poisson's ratio
parameters = {"MATRIX": [3.5, 0.35],   # Epoxy matrix
              "FIBER":  [230.0, 0.2]}  # Carbon fibers

# Name folder to store VTK files
vtk_path = os.path.join(workdir, "..", "vtk_output", str(fiber_fraction))

# ---------------------------------------------------------------
# Run homogenization routine
C = np.empty((3, 3))
for i in range(3):
    print("\n------------------------------")
    print("Run case {}.".format(i + 1))
    strain = np.array([0.0, 0.0, 0.0])
    strain[i] = 0.1 if i < 2 else 0.2
    print("Strain = " + str(strain.round(3).tolist()))

    # Create geometry in Gmsh and dump it to an input file
    model_name = "RVE_{}".format(i + 1)
    boundary_map = create_geometry(model_name, dx=dx, dy=dy,
                                   fiber_fraction=fiber_fraction)
    unit_cell = GmshRVE(model_name, boundary_map)
    unit_cell.generate_mesh(element_size=0.4)
    input_file = unit_cell.write_to_abaqus_input()

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
    strain_3d = np.array([strain[0], strain[1], 0.0, 0.0, 0.0, strain[2]])
    if boundary_conditions.upper() == "PERIODIC":
        unit_cell.apply_periodic_BC(step_name, strain_3d)
    elif boundary_conditions.upper() == "DIRICHLET":
        unit_cell.apply_affine_displacement_BC(step_name, strain_3d)
    else:
        raise RuntimeError("Invalid boundary conditions.")

    unit_cell.run_analysis(field_outputs=["NFORC", "U", "S"])
    stress = unit_cell.compute_averaged_stress(step_name)
    print("Stress = " + str(stress.round(3).tolist()))

    # Compute effective stiffness tensor
    C[:, i] = stress/strain[i]

    # Export to VTK
    unit_cell.export_to_vtk(vtk_path, field_outputs=["U", ("S", "MISES")],
                            sets=["ELSET_MATRIX", "ELSET_FIBER"])

gmsh.finalize()

print("\n------------------------------")
print("Effective constitutive matrix:")
print(C.round(3))

# Compute Reuss and Voigt bounds
C_matrix = isotropic_stiffness_tensor(*parameters["MATRIX"])
C_fiber = isotropic_stiffness_tensor(*parameters["FIBER"])
C_reuss = reuss(fiber_fraction, C_fiber, C_matrix)
C_voigt = voigt(fiber_fraction, C_fiber, C_matrix)

# Compute eigenvalues for comparison
eigen_reuss = np.sort(np.linalg.eig(C_reuss)[0])
eigen_voigt = np.sort(np.linalg.eig(C_voigt)[0])
eigen_numerical = np.sort(np.linalg.eig(C)[0])

# Export eigenvalues
file = str(fiber_fraction) + ".csv"
export_to_csv_file(eigen_reuss, os.path.join(result_path, "reuss", file))
export_to_csv_file(eigen_voigt, os.path.join(result_path, "voigt", file))
export_to_csv_file(eigen_voigt, os.path.join(result_path, "numerical", file))

print("\n------------------------------")
print("Eigenvalues:")
print("Voigt upper bound:     " + str(eigen_voigt.round(3).tolist()))
print("Numerical:             " + str(eigen_numerical.round(3).tolist()))
print("Reuss lower bound:     " + str(eigen_reuss.round(3).tolist()))
