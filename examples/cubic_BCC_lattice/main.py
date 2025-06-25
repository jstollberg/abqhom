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
from abqhom.abqmodel import AbaqusRVE, LoadHistory


def create_geometry(model_name, dx=1.0, dy=1.0, dz=1.0):
    # Initialize Gmsh and add a new model
    if not gmsh.is_initialized():
        gmsh.initialize()
    if model_name in gmsh.model.list():
        gmsh.model.set_current(model_name)
        gmsh.model.remove()
    gmsh.model.add(model_name)

    # Center point should be coordinate system origin
    x = y = z = 0.0

    # Add points on corners and in the center
    point_1 = gmsh.model.geo.add_point(x - 0.5*dx, y - 0.5*dy, z - 0.5*dz)
    point_2 = gmsh.model.geo.add_point(x + 0.5*dx, y - 0.5*dy, z - 0.5*dz)
    point_3 = gmsh.model.geo.add_point(x + 0.5*dx, y + 0.5*dy, z - 0.5*dz)
    point_4 = gmsh.model.geo.add_point(x - 0.5*dx, y + 0.5*dy, z - 0.5*dz)
    point_5 = gmsh.model.geo.add_point(x - 0.5*dx, y - 0.5*dy, z + 0.5*dz)
    point_6 = gmsh.model.geo.add_point(x + 0.5*dx, y - 0.5*dy, z + 0.5*dz)
    point_7 = gmsh.model.geo.add_point(x + 0.5*dx, y + 0.5*dy, z + 0.5*dz)
    point_8 = gmsh.model.geo.add_point(x - 0.5*dx, y + 0.5*dy, z + 0.5*dz)
    point_9 = gmsh.model.geo.add_point(x, y, z)

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

    # Add lines in the volume
    edge_13 = gmsh.model.geo.add_line(point_9, point_1)
    edge_14 = gmsh.model.geo.add_line(point_9, point_2)
    edge_15 = gmsh.model.geo.add_line(point_9, point_3)
    edge_16 = gmsh.model.geo.add_line(point_9, point_4)
    edge_17 = gmsh.model.geo.add_line(point_9, point_5)
    edge_18 = gmsh.model.geo.add_line(point_9, point_6)
    edge_19 = gmsh.model.geo.add_line(point_9, point_7)
    edge_20 = gmsh.model.geo.add_line(point_9, point_8)

    gmsh.model.geo.synchronize()

    # Add physical groups for primitive unit cells
    cubic_entities = [edge_1, edge_2, edge_3, edge_4, edge_5, edge_6,
                      edge_7, edge_8, edge_9, edge_10, edge_11, edge_12]
    cubic = gmsh.model.add_physical_group(1, cubic_entities)
    gmsh.model.set_physical_name(1, cubic, "CUBIC")
    bcc_entities = [edge_13, edge_14, edge_15, edge_16, edge_17, edge_18,
                    edge_19, edge_20]
    bcc = gmsh.model.add_physical_group(1, bcc_entities)
    gmsh.model.set_physical_name(1, bcc, "BCC")


# ---------------------------------------------------------------
# Set working directory
workdir = os.path.join(base_path,
                       "abqhom/examples/cubic_BCC_lattice/abq")
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
element_size = 1

# Material parameters: Young's modulus - Poisson's ratio - diameter
parameters = {"CUBIC": [210.0, 0.3, 0.5],
              "BCC": [210.0, 0.3, 1.0]}

# Multiply stiffness with weight factors due to overlaps on edges and faces
parameters["CUBIC"][0] *= 0.25

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
    model_name = "RUC_{}".format(i + 1)
    create_geometry(model_name, dx=dx, dy=dy, dz=dz)
    unit_cell = GmshRVE(model_name)
    unit_cell.generate_mesh(element_size=element_size)
    input_file = unit_cell.write_abaqus_input()

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
    unit_cell.apply_periodic_BC(step_name, history)
    unit_cell.run_analysis(field_outputs=["NFORC", "U"])
    stress = unit_cell.compute_averaged_stress(step_name)
    print("Stress = " + str(stress.round(4).tolist()))

    # Compute effective stiffness tensor
    C[:, i] = stress/strain_magnitude

    # Export to VTK
    unit_cell.export_to_vtk(
        vtk_path, field_outputs=["U"],
        radius={"ELSET_CUBIC": 0.5*parameters["CUBIC"][2],
                "ELSET_BCC": 0.5*parameters["BCC"][2]})

print("\n------------------------------")
print("Effective constitutive matrix:")
print(C.round(3))
