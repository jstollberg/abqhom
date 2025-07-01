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

    # Add points on surfaces and in center
    point_9 = gmsh.model.geo.add_point(x, y, z - 0.5*dz)
    point_10 = gmsh.model.geo.add_point(x, y, z + 0.5*dz)
    point_11 = gmsh.model.geo.add_point(x - 0.5*dx, y, z)
    point_12 = gmsh.model.geo.add_point(x + 0.5*dx, y, z)
    point_13 = gmsh.model.geo.add_point(x, y - 0.5*dy, z)
    point_14 = gmsh.model.geo.add_point(x, y + 0.5*dy, z)
    point_15 = gmsh.model.geo.add_point(x, y, z)

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

    # Add lines in volume
    edge_37 = gmsh.model.geo.add_line(point_15, point_1)
    edge_38 = gmsh.model.geo.add_line(point_15, point_2)
    edge_39 = gmsh.model.geo.add_line(point_15, point_3)
    edge_40 = gmsh.model.geo.add_line(point_15, point_4)
    edge_41 = gmsh.model.geo.add_line(point_15, point_5)
    edge_42 = gmsh.model.geo.add_line(point_15, point_6)
    edge_43 = gmsh.model.geo.add_line(point_15, point_7)
    edge_44 = gmsh.model.geo.add_line(point_15, point_8)

    gmsh.model.geo.synchronize()

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

    bcc_entities = [edge_37, edge_38, edge_39, edge_40, edge_41, edge_42,
                    edge_43, edge_44]
    bcc = gmsh.model.add_physical_group(1, bcc_entities)
    gmsh.model.set_physical_name(1, bcc, "BCC")


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
                       "abqhom/examples/test/abq")
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
dx = dy = dz = 15
element_size = 1.5

# Nonlinear geometry
nonlinear_geometry = True

# Boundary conditions
boundary_conditions = "PERIODIC"  # PERIODIC or DIRICHLET

# Material parameters+
parameters = {"CUBIC_ELASTIC":  [530e3, 0.49],
              "CUBIC_DIAMETER": [0.5],
              "FCC_ELASTIC":    [530e3, 0.49],
              "FCC_DIAMETER":   [0.5],
              "BCC_ELASTIC":    [530e3, 0.49],
              "BCC_DIAMETER":   [0.5],
              "DENSITY":        [1120]}

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
unit_cell.generate_mesh(element_order=2, element_size=element_size)

# Load the model into Abaqus
unit_cell = AbaqusRVE(unit_cell, custom_types={1: "B31"})

# Assign material
for i in ("CUBIC", "FCC", "BCC"):
    E, nu = parameters[i + "_ELASTIC"]
    rho = parameters["DENSITY"]

    material = unit_cell.model.Material(name=i)
    material.Elastic(table=[[E, nu]])
    material.Density(table=(rho, ))
    unit_cell.assign_beam_section(i, i, i,
                                  0.5*parameters[i + "_DIAMETER"][0])

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
    radius={i: 0.5*parameters[i + "_DIAMETER"][0] for i in ("CUBIC", "FCC",
                                                            "BCC")})

gmsh.finalize()
print("*** Finished ***")
