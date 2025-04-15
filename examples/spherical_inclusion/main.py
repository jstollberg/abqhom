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


def create_geometry(dx=1.0, dy=1.0, dz=1.0, radius=0.2):
    model_name = "RVE"

    # Initialize Gmsh and add a new model
    if not gmsh.is_initialized():
        gmsh.initialize()
    gmsh.model.add(model_name)

    gmsh.option.set_number("Geometry.Tolerance", 1e-5)

    # Center point should be coordinate system origin
    x = y = z = 0.0

    # Create the model
    gmsh.model.occ.add_box(x - dx/2, y - dy/2, z - dz/2, dx, dy, dz)
    gmsh.model.occ.add_sphere(x, y, z, radius)
    gmsh.model.occ.cut([(3, 1)], [(3, 2)])
    gmsh.model.occ.add_surface_loop([7])
    gmsh.model.occ.add_volume([1])
    gmsh.model.occ.synchronize()

    # Create map containing information on mesh periodicity and physical groups
    boundary_map = {"POINT_1": (0, [12]),
                    "POINT_2": (0, [15]),
                    "POINT_3": (0, [18]),
                    "POINT_4": (0, [14]),
                    "POINT_5": (0, [11]),
                    "POINT_6": (0, [16]),
                    "POINT_7": (0, [17]),
                    "POINT_8": (0, [13]),
                    "EDGE_1":  (1, [20]),
                    "EDGE_2":  (1, [27]),
                    "EDGE_3":  (1, [25]),
                    "EDGE_4":  (1, [19]),
                    "EDGE_5":  (1, [22]),
                    "EDGE_6":  (1, [24]),
                    "EDGE_7":  (1, [23]),
                    "EDGE_8":  (1, [17]),
                    "EDGE_9":  (1, [16]),
                    "EDGE_10": (1, [21]),
                    "EDGE_11": (1, [26]),
                    "EDGE_12": (1, [18]),
                    "FACE_1":  (2, [12]),
                    "FACE_2":  (2, [10]),
                    "FACE_3":  (2, [9]),
                    "FACE_4":  (2, [11]),
                    "FACE_5":  (2, [8]),
                    "FACE_6":  (2, [13])
                    }

    # Add physical groups
    matrix = gmsh.model.add_physical_group(3, [1])
    gmsh.model.set_physical_name(3, matrix, "MATRIX")
    matrix = gmsh.model.add_physical_group(3, [2])
    gmsh.model.set_physical_name(3, matrix, "INCLUSION")

    return model_name, boundary_map


# ---------------------------------------------------------------
# Set working directory
workdir = os.path.join(base_path,
                       "abqhom/examples/spherical_inclusion/abq")
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
dx = dy = dz = 10.0
radius = 0.4*dx
element_size = 1.2

# Material parameters: Young's modulus - Poisson's ratio
parameters = {"MATRIX": [10.0, 0.3],
              "INCLUSION":   [210.0, 0.3]}

# Create geometry in Gmsh and dump it to an input file
model_name, boundary_map = create_geometry(dx=dx, dy=dy, dz=dz, radius=radius)
unit_cell = GmshRVE(model_name, boundary_map)
unit_cell.generate_mesh(element_size=element_size)
input_file = unit_cell.write_to_abaqus_input()

# Load the model into Abaqus
unit_cell = AbaqusRVE(unit_cell, input_file)

# Assign material and beam profile sections
for name in parameters:
    E, nu = parameters[name]
    material = unit_cell.model.Material(name=name)
    material.Elastic(table=[[E, nu]])
    unit_cell.assign_solid_section(name, name, "ELSET_" + name)

# Generate assembly and step
step_name = "LOAD"
unit_cell.model.StaticStep(name=step_name, previous="Initial")

# Apply periodic boundary conditions and run analysis
strain = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
unit_cell.apply_periodic_BC(step_name, strain)
unit_cell.run_analysis(field_outputs=["NFORC", "U", "S", "E"])

# Export to VTK
unit_cell.export_to_vtk(vtk_path, field_outputs=["U", "S", "E"])
gmsh.finalize()
