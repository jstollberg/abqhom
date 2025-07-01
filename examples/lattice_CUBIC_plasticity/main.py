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
from abqhom.utils import voigt_to_tensor
from abqhom.gmshmodel import GmshRVE
from abqhom.abqmodel import AbaqusRVE, LoadHistory

# Import Abaqus stuff
from abaqusConstants import FIXED_USER_DEFINED_INC


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

    # Create geometry
    gmsh.model.occ.add_box(x - 0.5*dx, y - 0.5*dy, z - 0.5*dz, dx, dy, dz)
    gmsh.model.occ.synchronize()
    volumes, surfaces = gmsh.model.get_entities(3), gmsh.model.get_entities(2)
    gmsh.model.remove_entities(volumes)
    gmsh.model.remove_entities(surfaces)

    # Add physical group
    lines = gmsh.model.get_entities(1)
    cubic = gmsh.model.add_physical_group(1, [i for _, i in lines])
    gmsh.model.set_physical_name(1, cubic, "CUBIC")


def load_path(time):
    strain = np.zeros(6)
    strain[0] = np.sin(np.pi*time)
    strain[1] = np.exp(-((time - 0.5)/0.2)**2) - np.exp(-(0.5/0.2)**2)
    strain[2] = 4*(time - 0.5)**2 - 1
    strain[4] = np.sin(2*np.pi*time)**2
    strain[5] = np.sin(2*np.pi*time)
    strain *= 0.4
    strain = voigt_to_tensor(strain, "STRAIN")
    return strain


# ---------------------------------------------------------------
# Set working directory
workdir = os.path.join(base_path,
                       "abqhom/examples/lattice_CUBIC_plasticity/abq")
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
element_size = 10

# Nonlinear geometry
nonlinear_geometry = True

# Boundary conditions
boundary_conditions = "PERIODIC"  # PERIODIC or DIRICHLET

# Material parameters
parameters = {"ELASTIC":  [230.0, 0.2],
              "PLASTIC":  [250e-3, 5.0],
              "DENSITY":  [1e-5],
              "DIAMETER": [1.0]}

parameters["ELASTIC"][0] *= 0.25
parameters["PLASTIC"][0] *= 0.25
parameters["PLASTIC"][1] *= 0.25

# Create the load path for the strain
n_steps = 500
time = np.linspace(0.0, 1.0, n_steps)
history = LoadHistory(time, load_path)


# ---------------------------------------------------------------
# Create geometry in Gmsh and dump it to an input file
model_name = "RUC"
create_geometry(model_name, dx=dx, dy=dy, dz=dz)
unit_cell = GmshRVE(model_name)
unit_cell.add_periodic_mesh_constraints()
unit_cell.generate_mesh(element_order=1, element_size=element_size)
input_file = unit_cell.write_abaqus_input(custom_types={1: "B31"})

# Load the model into Abaqus
unit_cell = AbaqusRVE(unit_cell, custom_types={1: "B31"})

# Assign material
E, nu = parameters["ELASTIC"]
sig_y, H = parameters["PLASTIC"]
rho = parameters["DENSITY"]
material = unit_cell.model.Material(name="CUBIC")
material.Elastic(table=[[E, nu]])
material.Density(table=(rho, ))
material.Plastic(table=[[sig_y + H*eps_p, eps_p] for eps_p in
                        np.linspace(0.0, 5.0, 10)])
unit_cell.assign_beam_section("CUBIC", "CUBIC", "CUBIC",
                              0.5*parameters["DIAMETER"][0])

# Use an explicit scheme for stability reasons
step_name = "LOAD"
time_period = history.time_data[-1]
unit_cell.model.ExplicitDynamicsStep(
    name=step_name, previous="Initial",
    nlgeom=nonlinear_geometry,
    timePeriod=time_period,
    timeIncrementationMethod=FIXED_USER_DEFINED_INC,
    userDefinedInc=time_period/float(n_steps)
    )

# Apply boundary conditions and run analysis
if boundary_conditions.upper() == "PERIODIC":
    unit_cell.apply_periodic_BC(step_name, history)
elif boundary_conditions.upper() == "DIRICHLET":
    unit_cell.apply_affine_displacement_BC(step_name, history)
else:
    raise RuntimeError("Invalid boundary conditions.")

unit_cell.run_analysis(field_outputs=["NFORC", "U", "UR", "S"])

# Export to VTK
print("Running VTK export...")
unit_cell.export_to_vtk(
    vtk_path, field_outputs=["U", "UR"],
    radius={"CUBIC": 0.5*parameters["DIAMETER"][0]})

gmsh.finalize()
print("*** Finished ***")
