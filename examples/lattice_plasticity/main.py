"""
Elasto-plastic homogenization of a simple-cubic beam lattice under a
*cyclic* uniaxial load with growing amplitude (geometrically nonlinear,
static implicit).

This shows that the framework -- and the NFORC-based average_stress in
particular -- is not restricted to linear elasticity: the effective axial
response is evaluated frame by frame for an elasto-plastic strut material
with isotropic hardening. The macroscopic load is a strain-controlled cycle
whose amplitude grows a little every cycle, so the struts yield in tension
and compression and the (isotropically hardening) hysteresis loops nest
outward without crossing. The struts use their true material parameters;
redundant periodic copies of boundary struts are removed automatically
before the analysis.

Because the analysis is geometrically nonlinear (nlgeom=ON), the prescribed
macroscopic quantity is the displacement gradient H = F - I (F being the
macroscopic deformation gradient), not the linearized small-strain tensor.
Its axial component H_xx = F_xx - 1 (the nominal/engineering strain) is the
control variable, and it is work-conjugate to the first Piola-Kirchhoff
stress P_xx returned by average_stress. Both are plotted below.

The effective response is written to a CSV file and plotted with matplotlib
(bundled with Abaqus); the figure is saved next to the results.

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

from abaqusConstants import ON, ISOTROPIC


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

# Geometry and elasto-plastic material (bilinear isotropic hardening)
length = 10.0
radius = 0.5
element_size = 1.0
E, nu = 230.0, 0.3
yield_stress = 1.0        # strut material yields well below the peak load
hardening_modulus = 20.0  # post-yield tangent modulus

# Cyclic uniaxial macroscopic load, displacement-gradient controlled
# (H_xx = F_xx - 1). Fully reversed cycles whose amplitude grows every
# cycle, so the hysteresis loops nest outward instead of overlapping.
n_cycles = 3
H_max = 0.013                                    # amplitude of the last cycle
amplitudes = np.linspace(0.5*H_max, H_max, n_cycles)

# Reversal points of the path: 0, +A1, -A1, +A2, -A2, ... , 0
reversals = [0.0]
for a in amplitudes:
    reversals += [a, -a]
reversals.append(0.0)
reversals = np.array(reversals)
reversal_times = np.linspace(0.0, 1.0, len(reversals))

def load_path(time):
    H_xx = np.interp(time, reversal_times, reversals)
    gradient = np.zeros(6)
    gradient[0] = H_xx        # macroscopic H = F - I, axial component
    return voigt_to_tensor(gradient, "STRAIN")

history = LoadHistory(reversal_times, load_path)


# ---------------------------------------------------------------
# Create geometry and periodic mesh in Gmsh
model_name = "RUC"
create_geometry(model_name, length)
unit_cell = GmshRVE(model_name)
unit_cell.add_periodic_mesh_constraints()
unit_cell.generate_mesh(element_order=1, element_size=element_size)

# Load the model into Abaqus (2-node linear Timoshenko beams)
unit_cell = AbaqusRVE(unit_cell, custom_types={1: "B31"})

# Elasto-plastic material with isotropic hardening
material = unit_cell.model.Material(name="CUBIC")
material.Elastic(table=[[E, nu]])
material.Plastic(hardening=ISOTROPIC,
                 table=[[yield_stress + hardening_modulus*eps_p, eps_p]
                        for eps_p in np.linspace(0.0, 0.1, 5)])
unit_cell.assign_beam_section("CUBIC", "CUBIC", radius)

# Static, geometrically nonlinear step with automatic incrementation. The
# small maximum increment resolves the cyclic path with many output frames.
step_name = "LOAD"
unit_cell.model.StaticStep(name=step_name, previous="Initial", nlgeom=ON,
                           maxNumInc=1000, initialInc=0.01,
                           minInc=1.0e-6, maxInc=1.0/(20*n_cycles))

# Apply periodic boundary conditions and run the analysis
unit_cell.apply_periodic_BC(step_name, history)
unit_cell.run_analysis(field_outputs=["U", "S"])

# ---------------------------------------------------------------
# Effective axial response (H_xx, P_xx), evaluated frame by frame
unit_cell.open_odb_file()
frame_times = [frame.frameValue
               for frame in unit_cell.odb.steps[step_name].frames]
unit_cell.close_odb_file()

curve = []
for index, frame_time in enumerate(frame_times):
    # nlgeom=ON, so average_stress returns the first Piola-Kirchhoff stress
    # (not symmetric, hence not converted to Voigt notation); only its axial
    # component P_xx is needed. The macroscopic displacement gradient at this
    # frame, H_xx = F_xx - 1, is read back from the load path.
    H_xx = load_path(frame_time)[0, 0]
    stress = unit_cell.average_stress(step_name, frame=index)
    curve.append([H_xx, stress[0, 0]])
curve = np.array(curve)

# Peak tensile stress reached at each (growing) tension reversal -- it
# increases from cycle to cycle because of the isotropic hardening.
tension_peaks = reversal_times[1::2][:n_cycles]
print("\n------------------------------")
print("Cyclic effective axial response ({} cycles, {} frames):".format(
    n_cycles, len(curve)))
print("  initial (elastic) tangent = {:.4f}".format(curve[1, 1]/curve[1, 0]))
print("  cycle |  H_xx amplitude |  peak P_xx")
for c, (amp, t_peak) in enumerate(zip(amplitudes, tension_peaks)):
    frame = int(np.argmin(np.abs(np.array(frame_times) - t_peak)))
    print("  {:5d} |     {:.5f}     |  {:.6e}".format(c + 1, amp,
                                                      curve[frame, 1]))
print("  --> peak stress grows with the cycles: isotropic hardening.")

# Store the response curve and export the deformed lattice to VTK
result_dir = os.path.join(here, "results")
if not os.path.isdir(result_dir):
    os.makedirs(result_dir)
export_to_csv_file(curve, os.path.join(result_dir, "cyclic_response.csv"))

# ---------------------------------------------------------------
# Plot the cyclic stress/strain response (matplotlib ships with Abaqus).
# Use the non-interactive Agg backend so it works under 'noGUI'.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(6.0, 5.0))
ax.plot(curve[:, 0], curve[:, 1], "-", color="#1f6feb", lw=1.4)
ax.axhline(0.0, color="0.6", lw=0.6)
ax.axvline(0.0, color="0.6", lw=0.6)
ax.set_xlabel(r"macroscopic displacement gradient  $H_{xx} = F_{xx} - 1$")
ax.set_ylabel(r"effective first Piola-Kirchhoff stress  $P_{xx}$")
ax.set_title("Cyclic elasto-plastic response (isotropic hardening)")
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(result_dir, "cyclic_response.png"), dpi=150)
plt.close(fig)

unit_cell.export_to_vtk(os.path.join(here, "vtk_output"),
                        field_outputs=["U", ("S", "MISES")], clear=True)

gmsh.finalize()
print("*** Finished ***")
