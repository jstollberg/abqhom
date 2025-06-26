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
from abqhom.abqmodel import AbaqusRVE, LoadHistory
from abqhom.utils import reuss, voigt, export_to_csv_file

# Import Abaqus stuff
from abaqusConstants import FIXED


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

    # Construct basic geometry
    rectangle = gmsh.model.occ.add_rectangle(x - 0.5*dx, y - 0.5*dy, 0.0,
                                             dx, dy)
    circle_1 = gmsh.model.occ.add_disk(x, y, 0.0, radius, radius)
    circle_2 = gmsh.model.occ.add_disk(x - 0.5*dx, y - 0.5*dy, 0.0,
                                       radius, radius)
    circle_3 = gmsh.model.occ.add_disk(x + 0.5*dx, y - 0.5*dy, 0.0,
                                       radius, radius)
    circle_4 = gmsh.model.occ.add_disk(x + 0.5*dx, y + 0.5*dy, 0.0,
                                       radius, radius)
    circle_5 = gmsh.model.occ.add_disk(x - 0.5*dx, y + 0.5*dy, 0.0,
                                       radius, radius)

    # Perform boolean operations
    fibers, _ = gmsh.model.occ.intersect(
        [(2, circle_1), (2, circle_2), (2, circle_3), (2, circle_4),
         (2, circle_5)], [(2, rectangle)], removeTool=False)
    matrix, _ = gmsh.model.occ.cut([(2, rectangle)], fibers, removeTool=False)
    gmsh.model.occ.synchronize()

    # Add physical groups
    matrix = gmsh.model.add_physical_group(2, [i for _, i in matrix])
    gmsh.model.set_physical_name(2, matrix, "MATRIX")
    fibers = gmsh.model.add_physical_group(2, [i for _, i in fibers])
    gmsh.model.set_physical_name(2, fibers, "FIBER")


def load_path(time):
    strain = np.zeros(6)
    strain[0] = np.sin(np.pi*time)
    strain[1] = np.exp(-((time - 0.5)/0.2)**2) - np.exp(-(0.5/0.2)**2)
    strain[5] = np.sin(2*np.pi*time)
    strain *= 0.125
    return strain


# ---------------------------------------------------------------
# Set working directory
workdir = os.path.join(base_path,
                       "abqhom/examples/fiber_reinforced_2d_plasticity/abq")
if not os.path.isdir(workdir):
    os.mkdir(workdir)
os.chdir(workdir)

export_vtk = True

# RVE size
dx = dy = 10.0
fiber_fraction = 0.5  # Change the volume fraction of the fibers here

# Boundary conditions
boundary_conditions = "PERIODIC"  # PERIODIC or DIRICHLET

# Material parameters:
parameters = {"MATRIX_ELASTIC": [3.5, 0.35],
              "MATRIX_PLASTIC": [70e-3, 2.0],
              "FIBER_ELASTIC":  [230.0, 0.2],
              "FIBER_PLASTIC":  [250e-3, 5.0]}

# Name folder to store VTK files
vtk_path = os.path.join(workdir, "..", "vtk_output", str(fiber_fraction))

# Create the load path for the strain
n_steps = 500
time = np.linspace(0.0, 1.0, 50)
history = LoadHistory(time, load_path)


# ---------------------------------------------------------------
# Create geometry in Gmsh and dump it to an input file
model_name = "RVE"
create_geometry(model_name, dx=dx, dy=dy, fiber_fraction=fiber_fraction)
unit_cell = GmshRVE(model_name)
unit_cell.add_periodic_mesh_constraints()
unit_cell.generate_mesh(element_size=0.37)
input_file = unit_cell.write_abaqus_input()

# Load the model into Abaqus
unit_cell = AbaqusRVE(unit_cell, input_file)

# Assign material
for name in ("MATRIX", "FIBER"):
    E, nu = parameters[name + "_ELASTIC"]
    sig_y, H = parameters[name + "_PLASTIC"]
    material = unit_cell.model.Material(name=name)
    material.Elastic(table=[[E, nu]])
    material.Plastic(table=[[sig_y + H*eps_p, eps_p] for eps_p in
                            np.linspace(0.0, 5.0, 10)])
    unit_cell.assign_solid_section(name, name, name)

# Define load step
step_name = "LOAD"
time_period = history.time_data[-1]
unit_cell.model.StaticStep(name=step_name, previous="Initial",
                           timePeriod=time_period, maxNumInc=n_steps,
                           timeIncrementationMethod=FIXED,
                           initialInc=time_period/float(n_steps))

# Apply boundary conditions and run analysis
if boundary_conditions.upper() == "PERIODIC":
    unit_cell.apply_periodic_BC(step_name, history)
elif boundary_conditions.upper() == "DIRICHLET":
    unit_cell.apply_affine_displacement_BC(step_name, history)
else:
    raise RuntimeError("Invalid boundary conditions.")

unit_cell.run_analysis(field_outputs=["NFORC", "U", "S", "E", "PE", "PEEQ"])

# Export to VTK
if export_vtk:
    print("Running VTK export...")
    unit_cell.export_to_vtk(vtk_path,
                            field_outputs=["U", ("S", "MISES"), "PEEQ"],
                            element_sets=["MATRIX", "FIBER"])
print("*** Finished ***")
