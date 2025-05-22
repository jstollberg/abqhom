import os
import sys
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


def create_geometry(model_name, dx=1.0, dy=1.0, dz=1.0):
    # Initialize Gmsh and add a new model
    if not gmsh.is_initialized():
        gmsh.initialize()
    gmsh.model.add(model_name)

    # Center point should be coordinate system origin
    x = y = z = 0.0

    # Add points on corners
    point_1 = gmsh.model.geo.add_point(x - 0.5*dx, y - 0.5*dy, z - 0.5*dz)
    point_2 = gmsh.model.geo.add_point(x + 0.5*dx, y - 0.5*dy, z - 0.5*dz)
    point_3 = gmsh.model.geo.add_point(x + 0.5*dx, y + 0.5*dy, z - 0.5*dz)
    point_4 = gmsh.model.geo.add_point(x - 0.5*dx, y + 0.5*dy, z - 0.5*dz)
    point_5 = gmsh.model.geo.add_point(x - 0.5*dx, y - 0.5*dy, z + 0.5*dz)
    point_6 = gmsh.model.geo.add_point(x + 0.5*dx, y - 0.5*dy, z + 0.5*dz)
    point_7 = gmsh.model.geo.add_point(x + 0.5*dx, y + 0.5*dy, z + 0.5*dz)
    point_8 = gmsh.model.geo.add_point(x - 0.5*dx, y + 0.5*dy, z + 0.5*dz)

    # Add points on surfaces
    point_9 = gmsh.model.geo.add_point(x, y, z - 0.5*dz)
    point_10 = gmsh.model.geo.add_point(x, y, z + 0.5*dz)
    point_11 = gmsh.model.geo.add_point(x - 0.5*dx, y, z)
    point_12 = gmsh.model.geo.add_point(x + 0.5*dx, y, z)
    point_13 = gmsh.model.geo.add_point(x, y - 0.5*dy, z)
    point_14 = gmsh.model.geo.add_point(x, y + 0.5*dy, z)

    # Add lines on edges
    edge_1 = gmsh.model.geo.add_line(point_1, point_2)
    edge_2 = gmsh.model.geo.add_line(point_2, point_3)
    edge_3 = gmsh.model.geo.add_line(point_3, point_4)
    edge_4 = gmsh.model.geo.add_line(point_4, point_1)
    edge_5 = gmsh.model.geo.add_line(point_5, point_6)
    edge_6 = gmsh.model.geo.add_line(point_6, point_7)
    edge_7 = gmsh.model.geo.add_line(point_7, point_8)
    edge_8 = gmsh.model.geo.add_line(point_8, point_5)
    edge_9 = gmsh.model.geo.add_line(point_1, point_5)
    edge_10 = gmsh.model.geo.add_line(point_2, point_6)
    edge_11 = gmsh.model.geo.add_line(point_3, point_7)
    edge_12 = gmsh.model.geo.add_line(point_4, point_8)

    # Add lines on surfaces
    edge_13 = gmsh.model.geo.add_line(point_1, point_9)
    edge_14 = gmsh.model.geo.add_line(point_2, point_9)
    edge_15 = gmsh.model.geo.add_line(point_3, point_9)
    edge_16 = gmsh.model.geo.add_line(point_4, point_9)
    edge_17 = gmsh.model.geo.add_line(point_5, point_10)
    edge_18 = gmsh.model.geo.add_line(point_6, point_10)
    edge_19 = gmsh.model.geo.add_line(point_7, point_10)
    edge_20 = gmsh.model.geo.add_line(point_8, point_10)
    edge_21 = gmsh.model.geo.add_line(point_5, point_11)
    edge_22 = gmsh.model.geo.add_line(point_1, point_11)
    edge_23 = gmsh.model.geo.add_line(point_4, point_11)
    edge_24 = gmsh.model.geo.add_line(point_8, point_11)
    edge_25 = gmsh.model.geo.add_line(point_6, point_12)
    edge_26 = gmsh.model.geo.add_line(point_2, point_12)
    edge_27 = gmsh.model.geo.add_line(point_3, point_12)
    edge_28 = gmsh.model.geo.add_line(point_7, point_12)
    edge_29 = gmsh.model.geo.add_line(point_1, point_13)
    edge_30 = gmsh.model.geo.add_line(point_2, point_13)
    edge_31 = gmsh.model.geo.add_line(point_6, point_13)
    edge_32 = gmsh.model.geo.add_line(point_5, point_13)
    edge_33 = gmsh.model.geo.add_line(point_4, point_14)
    edge_34 = gmsh.model.geo.add_line(point_3, point_14)
    edge_35 = gmsh.model.geo.add_line(point_7, point_14)
    edge_36 = gmsh.model.geo.add_line(point_8, point_14)

    gmsh.model.geo.synchronize()

    # Create map containing information on mesh periodicity and physical groups
    boundary_map = {"POINT_1": (0, [point_1]),
                    "POINT_2": (0, [point_2]),
                    "POINT_3": (0, [point_3]),
                    "POINT_4": (0, [point_4]),
                    "POINT_5": (0, [point_5]),
                    "POINT_6": (0, [point_6]),
                    "POINT_7": (0, [point_7]),
                    "POINT_8": (0, [point_8]),
                    "EDGE_1": (1, [edge_1]),
                    "EDGE_2": (1, [edge_2]),
                    "EDGE_3": (1, [edge_3]),
                    "EDGE_4": (1, [edge_4]),
                    "EDGE_5": (1, [edge_5]),
                    "EDGE_6": (1, [edge_6]),
                    "EDGE_7": (1, [edge_7]),
                    "EDGE_8": (1, [edge_8]),
                    "EDGE_9": (1, [edge_9]),
                    "EDGE_10": (1, [edge_10]),
                    "EDGE_11": (1, [edge_11]),
                    "EDGE_12": (1, [edge_12]),
                    "FACE_1": (1, [edge_13, edge_14, edge_15, edge_16]),
                    "FACE_2": (1, [edge_17, edge_18, edge_19, edge_20]),
                    "FACE_3": (1, [edge_29, edge_30, edge_31, edge_32]),
                    "FACE_4": (1, [edge_33, edge_34, edge_35, edge_36]),
                    "FACE_5": (1, [edge_21, edge_22, edge_23, edge_24]),
                    "FACE_6": (1, [edge_25, edge_26, edge_27, edge_28])
                    }

    # Add physical groups for primitive unit cells
    cubic_entities = [edge_1, edge_2, edge_3, edge_4, edge_5, edge_6,
                      edge_7, edge_8, edge_9, edge_10, edge_11, edge_12]
    cubic = gmsh.model.add_physical_group(1, cubic_entities)
    gmsh.model.set_physical_name(1, cubic, "CUBIC")

    fcc_entities = [edge_13, edge_14, edge_15, edge_16, edge_17, edge_18,
                    edge_19, edge_20, edge_21, edge_22, edge_23, edge_24,
                    edge_25, edge_26, edge_27, edge_28, edge_29, edge_30,
                    edge_31, edge_32, edge_33, edge_34, edge_35, edge_36]
    fcc = gmsh.model.add_physical_group(1, fcc_entities)
    gmsh.model.set_physical_name(1, fcc, "FCC")

    return boundary_map


# ---------------------------------------------------------------
# Set working directory
workdir = os.path.join(base_path,
                       "abqhom/examples/cubic_FCC_lattice/abq")
if not os.path.isdir(workdir):
    os.mkdir(workdir)
os.chdir(workdir)

# Create folder to store results
result_path = os.path.join(workdir, "..", "results")
if not os.path.isdir(result_path):
    os.mkdir(result_path)

# Name folder to store VTK files
vtk_path = os.path.join(workdir, "..", "vtk_output")

# RVE size
dx = dy = dz = 10
element_size = 0.5

# Material parameters: Young's modulus - Poisson's ratio - diameter
parameters = {"CUBIC": [210.0, 0.3, 0.3],
              "FCC":   [210.0, 0.3, 0.6]}

# Multiply stiffness with weight factors due to overlaps on edges and faces
parameters["CUBIC"][0] *= 0.25
parameters["FCC"][0] *= 0.5

# Run homogenization routine
C = np.empty((6, 6))
for i in range(6):
    print("\n------------------------------")
    print("Run case {}.".format(i + 1))
    strain = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    strain[i] = 0.1 if i < 3 else 0.2
    print("Strain = " + str(strain.round(3).tolist()))

    # Create geometry in Gmsh and dump it to an input file
    model_name = "RUC_{}".format(i + 1)
    boundary_map = create_geometry(model_name, dx=dx, dy=dy, dz=dz)
    unit_cell = GmshRVE(model_name, boundary_map)
    unit_cell.generate_mesh(element_size=element_size)
    input_file = unit_cell.write_to_abaqus_input()

    # Load the model into Abaqus
    unit_cell = AbaqusRVE(unit_cell, input_file)
    gmsh.finalize()

    # Assign material
    for name in parameters:
        E, nu, d = parameters[name]
        material = unit_cell.model.Material(name=name)
        material.Elastic(table=[[E, nu]])
        unit_cell.assign_beam_section(name, name, "ELSET_" + name, 0.5*d)

    # Define load step
    step_name = "LOAD"
    unit_cell.model.StaticStep(name=step_name, previous="Initial")

    # Apply periodic boundary conditions and run analysis
    unit_cell.apply_periodic_BC(step_name, strain)
    unit_cell.run_analysis(field_outputs=["NFORC", "U"])
    stress = unit_cell.compute_averaged_stress(step_name)
    print("Stress = " + str(stress.round(4).tolist()))

    # Compute effective stiffness tensor
    C[:, i] = stress/strain[i]

    # Export to VTK
    unit_cell.export_to_vtk(
        vtk_path, field_outputs=["U"],
        radius={"ELSET_CUBIC": 0.5*parameters["CUBIC"][2],
                "ELSET_FCC": 0.5*parameters["FCC"][2]})

print("\n------------------------------")
print("Effective constitutive matrix:")
print(C.round(3))

print >> sys.__stdout__, "\n------------------------------"
print >> sys.__stdout__, "Effective constitutive matrix:"
print >> sys.__stdout__, C.round(3)
