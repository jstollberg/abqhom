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
from abqhom.utils import voigt_to_tensor, tensor_to_voigt
from abqhom.gmshmodel import GmshRVE
from abqhom.abqmodel import AbaqusRVE, LoadHistory
from abqhom.utils import reuss, voigt, export_to_csv_file


def create_geometry(model_name, dx=1.0, dy=1.0, dz=1.0, fiber_fraction=0.2):
    # Compute fiber radius from volume fraction
    radius = np.sqrt((fiber_fraction*dx*dy)/(2*np.pi))
    if radius >= np.sqrt((dx/2.0)**2.0 + (dy/2.0)**2.0)/2.0:
        raise ValueError("Fibers overlap.")

    # Initialize Gmsh and add a new model
    if not gmsh.is_initialized():
        gmsh.initialize()
    gmsh.model.add(model_name)

    gmsh.option.set_number("Geometry.Tolerance", 1e-5)

    # Center point should be coordinate system origin
    x = y = z = 0.0

    # Construct basic geometry
    rectangle = gmsh.model.occ.add_rectangle(x - 0.5*dx, y - 0.5*dy,
                                             z - 0.5*dz, dx, dy)
    circle_1 = gmsh.model.occ.add_disk(x, y, z - 0.5*dz, radius, radius)
    circle_2 = gmsh.model.occ.add_disk(x - 0.5*dx, y - 0.5*dy, z - 0.5*dz,
                                       radius, radius)
    circle_3 = gmsh.model.occ.add_disk(x + 0.5*dx, y - 0.5*dy, z - 0.5*dz,
                                       radius, radius)
    circle_4 = gmsh.model.occ.add_disk(x + 0.5*dx, y + 0.5*dy, z - 0.5*dz,
                                       radius, radius)
    circle_5 = gmsh.model.occ.add_disk(x - 0.5*dx, y + 0.5*dy, z - 0.5*dz,
                                       radius, radius)

    # Perform boolean operations
    fibers, _ = gmsh.model.occ.intersect(
        [(2, circle_1), (2, circle_2), (2, circle_3), (2, circle_4),
         (2, circle_5)], [(2, rectangle)], removeTool=False)
    matrix, _ = gmsh.model.occ.cut([(2, rectangle)], fibers, removeTool=False)

    # Extrude in z-direction
    fibers = gmsh.model.occ.extrude(fibers, 0.0, 0.0, dz)
    matrix = gmsh.model.occ.extrude(matrix, 0.0, 0.0, dz)
    gmsh.model.occ.synchronize()

    # Add physical groups
    matrix = gmsh.model.add_physical_group(3, [i for _, i in matrix])
    gmsh.model.set_physical_name(3, matrix, "MATRIX")
    fibers = gmsh.model.add_physical_group(3, [i for _, i in fibers])
    gmsh.model.set_physical_name(3, fibers, "FIBER")


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
                       "abqhom/examples/fiber_reinforced_3d/abq")
if not os.path.isdir(workdir):
    os.mkdir(workdir)
os.chdir(workdir)

debug_geometry = False

# RVE size
dx = dy = dz = 10.0
min_fiber_fraction = 0.5
max_fiber_fraction = 1.0
samples = 25

# Boundary conditions
boundary_conditions = "DIRICHLET"  # PERIODIC or DIRICHLET

# Material parameters for CFRP: Young's modulus - Poisson's ratio
parameters = {"MATRIX": [3.5, 0.35],   # Epoxy matrix
              "FIBER":  [230.0, 0.2]}  # Carbon fibers

# Create folders to store results
result_path = os.path.join(workdir, "..", "results", boundary_conditions)
if not os.path.isdir(result_path):
    os.makedirs(os.path.join(result_path, "reuss"))
    os.makedirs(os.path.join(result_path, "voigt"))
    os.makedirs(os.path.join(result_path, "numerical"))

# Test for different fiber volume fractions
for fiber_fraction in np.linspace(min_fiber_fraction,
                                  max_fiber_fraction,
                                  samples):
    print("\n\n================================")
    print("Fiber volume fraction = " + str(fiber_fraction) + "\n")

    # Name folder to store VTK files
    vtk_path = os.path.join(workdir, "..", "vtk_output", boundary_conditions,
                            str(fiber_fraction))

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
        create_geometry(model_name, dx=dx, dy=dy, dz=dz,
                        fiber_fraction=fiber_fraction)
        unit_cell = GmshRVE(model_name)
        unit_cell.add_periodic_mesh_constraints()
        unit_cell.generate_mesh(element_size=1.25)
        input_file = unit_cell.write_abaqus_input()

        if debug_geometry:
            unit_cell.show()
            sys.exit()

        # Load the model into Abaqus
        unit_cell = AbaqusRVE(unit_cell, input_file)
        gmsh.finalize()

        # Assign material
        for name in parameters:
            E, nu = parameters[name]
            material = unit_cell.model.Material(name=name)
            material.Elastic(table=[[E, nu]])
            unit_cell.assign_solid_section(name, name, name)

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
        stress = tensor_to_voigt(stress, "STRESS")
        print("Stress = " + str(stress.round(3).tolist()))

        # Compute effective stiffness tensor
        C[:, i] = stress/strain_magnitude

        # Export to VTK
        unit_cell.export_to_vtk(vtk_path, field_outputs=["U", ("S", "MISES")],
                                element_sets=["MATRIX", "FIBER"])

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
    export_to_csv_file(eigen_numerical, os.path.join(result_path, "numerical",
                                                     file))

    print("\n------------------------------")
    print("Eigenvalues:")
    print("Voigt upper bound:     " + str(eigen_voigt.round(3).tolist()))
    print("Numerical:             " + str(eigen_numerical.round(3).tolist()))
    print("Reuss lower bound:     " + str(eigen_reuss.round(3).tolist()))
