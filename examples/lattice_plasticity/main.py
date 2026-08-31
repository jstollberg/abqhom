"""
Elasto-plastic simulation of a simple-cubic beam lattice under a cyclic
multiaxial strain path (geometrically nonlinear, explicit dynamics).

This example shows that the framework is not limited to linear elasticity:
an arbitrary nonlinear material and load history can be prescribed. Note
that the struts use their true material parameters -- no manual stiffness
scaling is needed, because periodic strut sharing is handled downstream in
average_stress (which is material-agnostic) rather than by modifying the
material here.

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

from abqhom.utils import voigt_to_tensor
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
# Working directory (created next to this script)
workdir = os.path.join(here, "abq")
if not os.path.isdir(workdir):
    os.makedirs(workdir)
os.chdir(workdir)

vtk_path = os.path.join(here, "vtk_output")

# RVE size and mesh
dx = dy = dz = 10
element_size = 10
nonlinear_geometry = True

# Material parameters (true strut material -- no sharing scaling here)
parameters = {"ELASTIC":  [230.0, 0.2],   # Young's modulus - Poisson's ratio
              "PLASTIC":  [250e-3, 5.0],  # yield stress - hardening modulus
              "DENSITY":  [1e-5],
              "DIAMETER": [1.0]}

# Create the cyclic strain load path
n_steps = 500
time = np.linspace(0.0, 1.0, n_steps)
history = LoadHistory(time, load_path)


# ---------------------------------------------------------------
# Create geometry and periodic mesh in Gmsh
model_name = "RUC"
create_geometry(model_name, dx=dx, dy=dy, dz=dz)
unit_cell = GmshRVE(model_name)
unit_cell.add_periodic_mesh_constraints()
unit_cell.generate_mesh(element_order=1, element_size=element_size)

# Load the model into Abaqus (2-node linear Timoshenko beams)
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
unit_cell.assign_beam_section("CUBIC", "CUBIC", 0.5*parameters["DIAMETER"][0])

# Use an explicit scheme for robustness. The stable time increment is
# chosen automatically by Abaqus from the actual material stiffness, so no
# manual increment (or material softening) is required.
step_name = "LOAD"
time_period = history.time_data[-1]
unit_cell.model.ExplicitDynamicsStep(
    name=step_name, previous="Initial",
    nlgeom=nonlinear_geometry,
    timePeriod=time_period)

# Apply periodic boundary conditions and run the analysis
unit_cell.apply_periodic_BC(step_name, history)
unit_cell.run_analysis(field_outputs=["U", "UR", "S"])

# Export to VTK (the strut radius is taken from the assigned section)
print("Running VTK export...")
unit_cell.export_to_vtk(vtk_path, field_outputs=["U", "UR"], clear=True)

gmsh.finalize()
print("*** Finished ***")
