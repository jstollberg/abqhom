"""
Elasto-plastic homogenization of a simple-cubic beam lattice under a
monotonic uniaxial strain (geometrically nonlinear, static implicit).

This shows that the framework -- and the NFORC-based average_stress in
particular -- is not restricted to linear elasticity: the effective axial
stress/strain response is evaluated frame by frame for an elasto-plastic
strut material, exhibiting the transition from the elastic slope to
plastic yielding. The struts use their true material parameters; periodic
strut sharing is handled automatically in average_stress.

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

from abqhom.utils import voigt_to_tensor, export_to_csv_file
from abqhom.gmshmodel import GmshRVE
from abqhom.abqmodel import AbaqusRVE, LoadHistory

from abaqusConstants import ON


def create_geometry(model_name, length):
    """Create the twelve edge struts of a cube of the given edge length."""
    if not gmsh.is_initialized():
        gmsh.initialize()
    if model_name in gmsh.model.list():
        gmsh.model.set_current(model_name)
        gmsh.model.remove()
    gmsh.model.add(model_name)
    g = gmsh.model.geo

    h = 0.5*length
    coords = [(x, y, z) for x in (-h, h) for y in (-h, h) for z in (-h, h)]
    corners = [g.add_point(*c) for c in coords]

    struts = []
    for i in range(8):
        for j in range(i + 1, 8):
            shared = 0
            for a, b in zip(coords[i], coords[j]):
                if a == b:
                    shared += 1
            if shared == 2:
                struts.append(g.add_line(corners[i], corners[j]))

    g.synchronize()
    group = gmsh.model.add_physical_group(1, struts)
    gmsh.model.set_physical_name(1, group, "CUBIC")


# ---------------------------------------------------------------
# Working directory (created next to this script)
workdir = os.path.join(here, "abq")
if not os.path.isdir(workdir):
    os.makedirs(workdir)
os.chdir(workdir)

# Geometry and elasto-plastic material (bilinear hardening)
length = 10.0
radius = 0.5
element_size = 1.0
E, nu = 230.0, 0.3
yield_stress = 1.0        # strut material yields well below the peak load
hardening_modulus = 20.0  # post-yield tangent modulus

# Monotonic uniaxial macroscopic strain, ramped from 0 to eps_max
eps_max = 0.01

def load_path(time):
    strain = np.zeros(6)
    strain[0] = eps_max*time
    return voigt_to_tensor(strain, "STRAIN")

history = LoadHistory(np.linspace(0.0, 1.0, 25), load_path)


# ---------------------------------------------------------------
# Create geometry and periodic mesh in Gmsh
model_name = "RUC"
create_geometry(model_name, length)
unit_cell = GmshRVE(model_name)
unit_cell.add_periodic_mesh_constraints()
unit_cell.generate_mesh(element_order=1, element_size=element_size)

# Load the model into Abaqus (2-node linear Timoshenko beams)
unit_cell = AbaqusRVE(unit_cell, custom_types={1: "B31"})

# Elasto-plastic material
material = unit_cell.model.Material(name="CUBIC")
material.Elastic(table=[[E, nu]])
material.Plastic(table=[[yield_stress + hardening_modulus*eps_p, eps_p]
                        for eps_p in np.linspace(0.0, 0.1, 5)])
unit_cell.assign_beam_section("CUBIC", "CUBIC", radius)

# Static, geometrically nonlinear step with automatic incrementation. The
# small maximum increment produces several output frames along the path.
step_name = "LOAD"
unit_cell.model.StaticStep(name=step_name, previous="Initial", nlgeom=ON,
                           maxNumInc=200, initialInc=0.05, minInc=1.0e-6,
                           maxInc=0.05)

# Apply periodic boundary conditions and run the analysis
unit_cell.apply_periodic_BC(step_name, history)
unit_cell.run_analysis(field_outputs=["U", "S"])

# ---------------------------------------------------------------
# Effective axial stress/strain response, evaluated frame by frame
unit_cell.open_odb_file()
frame_times = [frame.frameValue
               for frame in unit_cell.odb.steps[step_name].frames]
unit_cell.close_odb_file()

curve = []
for index, frame_time in enumerate(frame_times):
    # The nonlinear analysis yields the first Piola-Kirchhoff stress, which
    # is not symmetric and hence not converted to Voigt notation; only the
    # axial component P_xx is needed for the effective response.
    stress = unit_cell.average_stress(step_name, frame=index)
    curve.append([eps_max*frame_time, stress[0, 0]])
curve = np.array(curve)

print("\n------------------------------")
print("Effective axial response (strain_xx, stress_xx):")
for eps, sig in curve:
    print("  {:.5f}   {:.6e}".format(eps, sig))
print("\ninitial (elastic) tangent C11 = {:.4f}".format(
    curve[1, 1]/curve[1, 0]))
print("final   (plastic) tangent     = {:.4f}".format(
    (curve[-1, 1] - curve[-2, 1])/(curve[-1, 0] - curve[-2, 0])))

# Store the response curve and export the deformed lattice to VTK
result_dir = os.path.join(here, "results")
if not os.path.isdir(result_dir):
    os.makedirs(result_dir)
export_to_csv_file(curve, os.path.join(result_dir, "axial_response.csv"))

unit_cell.export_to_vtk(os.path.join(here, "vtk_output"),
                        field_outputs=["U", ("S", "MISES")], clear=True)

gmsh.finalize()
print("*** Finished ***")
