"""
Linear homogenization of a 2d regular honeycomb beam-lattice unit cell.

This example demonstrates that abqhom now handles general lattice unit
cells whose struts do NOT connect to the corners of the bounding box.
The regular honeycomb is a classic case: its cell walls cross the
bounding-box edges at non-corner, periodically matching points.

Note that the workflow is *identical* to the one for classic "cubic"
cells -- no corner points and no manual boundary or rigid-body-anchor
definitions are required. The boundary, its periodic connections and the
node used to suppress rigid-body motion are all detected automatically.

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

from abqhom.utils import voigt_to_tensor, tensor_to_voigt, export_to_csv_file
from abqhom.gmshmodel import GmshRVE
from abqhom.abqmodel import AbaqusRVE, LoadHistory


def _clip_segment(p, q, box, tol=1e-9):
    """Liang-Barsky clip of the segment p->q to box=(xmin,xmax,ymin,ymax)."""
    xmin, xmax, ymin, ymax = box
    x0, y0 = p
    x1, y1 = q
    dx, dy = x1 - x0, y1 - y0
    t0, t1 = 0.0, 1.0
    for pk, qk in ((-dx, x0 - xmin), (dx, xmax - x0),
                   (-dy, y0 - ymin), (dy, ymax - y0)):
        if abs(pk) < tol:
            if qk < -tol:
                return None
        else:
            r = qk/pk
            if pk < 0:
                if r > t1:
                    return None
                t0 = max(t0, r)
            else:
                if r < t0:
                    return None
                t1 = min(t1, r)
    if t1 - t0 < tol:
        return None
    return ((x0 + t0*dx, y0 + t0*dy), (x0 + t1*dx, y0 + t1*dy))


def create_honeycomb(model_name, a=1.0, offset=(0.31, 0.17)):
    """Create a regular honeycomb beam lattice in a rectangular periodic
    unit cell (3a x sqrt(3)a) by clipping the infinite honeycomb network
    to the box. Clipping is translation-covariant, so opposite boundaries
    match and the cell is periodic without any node on a box corner.
    """
    if not gmsh.is_initialized():
        gmsh.initialize()
    if model_name in gmsh.model.list():
        gmsh.model.set_current(model_name)
        gmsh.model.remove()
    gmsh.model.add(model_name)
    g = gmsh.model.geo

    s3 = np.sqrt(3.0)
    Lx, Ly = 3.0*a, s3*a
    ox, oy = offset
    box = (ox, ox + Lx, oy, oy + Ly)

    a1 = np.array([1.5, s3/2.0])*a
    a2 = np.array([1.5, -s3/2.0])*a
    deltas = [np.array([1.0, 0.0])*a,
              np.array([-0.5, s3/2.0])*a,
              np.array([-0.5, -s3/2.0])*a]

    # Collect clipped strut segments over a window covering the box
    segments = []
    for i in range(-4, 5):
        for j in range(-4, 5):
            A = i*a1 + j*a2
            for d in deltas:
                clipped = _clip_segment(tuple(A), tuple(A + d), box)
                if clipped is not None:
                    segments.append(clipped)

    # Deduplicate points and lines
    def key(pt):
        return (round(pt[0], 7), round(pt[1], 7))

    point_tags = {}

    def get_point(pt):
        k = key(pt)
        if k not in point_tags:
            point_tags[k] = g.add_point(k[0], k[1], 0.0)
        return point_tags[k]

    seen = set()
    lines = []
    for p, q in segments:
        kp, kq = key(p), key(q)
        if kp == kq:
            continue
        ek = tuple(sorted((kp, kq)))
        if ek in seen:
            continue
        seen.add(ek)
        lines.append(g.add_line(get_point(p), get_point(q)))

    g.synchronize()
    grp = gmsh.model.add_physical_group(1, lines)
    gmsh.model.set_physical_name(1, grp, "STRUTS")


# ---------------------------------------------------------------
# Working directory (created next to this script)
workdir = os.path.join(here, "abq")
if not os.path.isdir(workdir):
    os.makedirs(workdir)
os.chdir(workdir)

# Material parameters (of the solid strut material)
E = 1.0     # Young's modulus
nu = 0.3    # Poisson's ratio

# Lattice parameters
a = 1.0             # cell-wall length
radius = 0.05       # strut cross-section radius
element_size = 0.2  # beam element length

# ---------------------------------------------------------------
# Run the homogenization routine (3 load cases in 2d: xx, yy, xy)
components = [0, 1, 5]  # Voigt indices: xx, yy, xy(engineering)
strain_magnitude = 1.0e-3

C = np.empty((3, 3))
for col, comp in enumerate(components):
    print("\n------------------------------")
    print("Run load case {}.".format(col + 1))

    def strain_history(time):
        strain = np.zeros(6)
        strain[comp] = strain_magnitude
        return time*voigt_to_tensor(strain, "STRAIN")

    history = LoadHistory(np.array([0.0, 1.0]), strain_history)

    # Create geometry and periodic mesh in Gmsh
    model_name = "RVE_{}".format(col + 1)
    create_honeycomb(model_name, a=a)
    unit_cell = GmshRVE(model_name)
    unit_cell.add_periodic_mesh_constraints()
    unit_cell.generate_mesh(element_size=element_size)

    # Load the model into Abaqus and assign a beam section to the struts
    unit_cell = AbaqusRVE(unit_cell)
    material = unit_cell.model.Material(name="LINEAR_ELASTICITY")
    material.Elastic(table=[[E, nu]])
    unit_cell.assign_beam_section("STRUTS", "LINEAR_ELASTICITY", radius)

    # Add load step and periodic boundary conditions (anchor is picked
    # automatically -- no corner point required)
    unit_cell.model.StaticStep(name="LOAD", previous="Initial")
    unit_cell.apply_periodic_BC("LOAD", history)

    # Run the analysis and compute the average stress
    unit_cell.run_analysis(field_outputs=["U"])
    stress = unit_cell.average_stress("LOAD")
    stress = tensor_to_voigt(stress, "STRESS", check_symmetry=False)
    print("Stress = " + str(stress.round(8).tolist()))

    C[:, col] = stress/strain_magnitude

    # Export the deformed struts to VTK (the strut radius is taken
    # automatically from the assigned beam section)
    unit_cell.export_to_vtk(os.path.join(here, "vtk_output"),
                            field_outputs=["U"], clear=(col == 0))
    gmsh.finalize()

print("\n------------------------------")
print("Effective constitutive matrix (Voigt):")
print(C)
print("\nChecks for the regular honeycomb (should be in-plane isotropic):")
print("  symmetry   max|C - C^T| = {:.2e}".format(np.max(np.abs(C - C.T))))
print("  C11 - C22               = {:.2e}".format(C[0, 0] - C[1, 1]))
print("  (C11 - C12) - 2 C33     = {:.2e}".format(
    (C[0, 0] - C[0, 1]) - 2.0*C[2, 2]))

# Export the effective stiffness
result_path = os.path.join(here, "results")
if not os.path.isdir(result_path):
    os.makedirs(result_path)
export_to_csv_file(C, os.path.join(result_path, "honeycomb_C.csv"))
