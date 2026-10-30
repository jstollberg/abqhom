"""Homogenize a re-entrant shell honeycomb in the linear regime.

The ribs of a two-dimensional re-entrant honeycomb in the x--z plane are
extruded along y to form thin shell walls. The infinite honeycomb is clipped
to one periodic cell before extrusion. The example computes the complete
effective stiffness and reports the in-plane auxetic response.
"""

import inspect
import os
import sys

import numpy as np

# ---------------------------------------------------------------
# Locate this example and put the shipped abqhom source first on the
# path. Under 'abaqus cae noGUI=main.py' __file__ is undefined, so the
# script location is obtained via inspect instead.
here = os.path.dirname(os.path.abspath(
    inspect.getframeinfo(inspect.currentframe()).filename))
abqhom_path = os.path.abspath(os.path.join(here, os.pardir, os.pardir, "src"))
if abqhom_path not in sys.path:
    sys.path.insert(0, abqhom_path)

import gmsh

from abqhom.abqmodel import AbaqusRVE, LoadHistory
from abqhom.constants import PARAQUS_INSTALLED
from abqhom.gmshmodel import GmshRVE
from abqhom.utils import export_to_csv_file, tensor_to_voigt, voigt_to_tensor


def _clip_segment(a, b, box, tol=1e-9):
    """Clip an in-plane segment to a box with the Liang--Barsky method."""
    x0, z0 = a
    x1, z1 = b
    dx, dz = x1 - x0, z1 - z0
    t0, t1 = 0.0, 1.0
    for pk, qk in ((-dx, x0 - box[0]), (dx, box[1] - x0),
                   (-dz, z0 - box[2]), (dz, box[3] - z0)):
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
    return ((x0 + t0*dx, z0 + t0*dz), (x0 + t1*dx, z0 + t1*dz))


def create_reentrant_honeycomb(model_name, h=1.0, l=1.0, theta_deg=-18.0,
                               depth=1.0, offset=(0.37, 0.29)):
    """Create a periodic re-entrant honeycomb from shell walls.

    The 2D honeycomb ribs (vertical ribs of length ``h`` plus inclined ribs
    of length ``l`` at the re-entrant angle ``theta_deg`` < 0) are laid out
    on the node set {(i, k) : i + k even}, clipped to one in-plane period
    (2*l*cos(theta) x 2*(h + l*sin(theta))), and each rib is extruded a
    depth ``depth`` along y into a plate wall (physical group "WALLS").
    """
    if theta_deg >= 0.0:
        raise ValueError("Need theta < 0 for a re-entrant (auxetic) cell.")

    theta = np.radians(theta_deg)
    p = l*np.cos(theta)
    q = h + l*np.sin(theta)            # in-plane half-period (< h)

    if not gmsh.is_initialized():
        gmsh.initialize()
    if model_name in gmsh.model.list():
        gmsh.model.set_current(model_name)
        gmsh.model.remove()
    gmsh.model.add(model_name)
    g = gmsh.model.geo

    ox, oz = offset[0]*p, offset[1]*q
    box = (ox, ox + 2*p, oz, oz + 2*q)

    # Rib segments of the infinite honeycomb, clipped to one period
    segments = []
    for i in range(-4, 6):
        for k in range(-4, 8):
            if (i + k) % 2 != 0:
                continue
            x, z = i*p, k*q
            for a, b in (((x, z), (x, z + h)),               # vertical rib
                         ((x, z + h), (x + p, z + q)),       # inclined ribs
                         ((x, z + h), (x - p, z + q))):
                clipped = _clip_segment(a, b, box)
                if clipped is not None:
                    segments.append(clipped)

    # Deduplicated points and (signed) lines; each rib is extruded into a
    # rectangular plate spanning y in [0, depth].
    point_tags = {}

    def get_point(x, y, z):
        """Return the coordinate key of a deduplicated geometric point."""
        key = (round(x, 6), round(y, 6), round(z, 6))
        if key not in point_tags:
            point_tags[key] = g.add_point(*key)
        return key

    line_tags = {}

    def get_line(ka, kb):
        """Return an oriented line tag between two point keys."""
        key = tuple(sorted((ka, kb)))
        if key not in line_tags:
            line_tags[key] = g.add_line(point_tags[key[0]],
                                        point_tags[key[1]])
        return line_tags[key] if key[0] == ka else -line_tags[key]

    surfaces = []
    seen = set()
    for (x0, z0), (x1, z1) in segments:
        ka2 = (round(x0, 6), round(z0, 6))
        kb2 = (round(x1, 6), round(z1, 6))
        if ka2 == kb2:
            continue
        ek = tuple(sorted((ka2, kb2)))
        if ek in seen:
            continue
        seen.add(ek)
        a0 = get_point(x0, 0.0, z0)
        a1 = get_point(x1, 0.0, z1)
        b1 = get_point(x1, depth, z1)
        b0 = get_point(x0, depth, z0)
        loop = g.add_curve_loop([get_line(a0, a1), get_line(a1, b1),
                                 get_line(b1, b0), get_line(b0, a0)])
        surfaces.append(g.add_plane_surface([loop]))

    g.synchronize()
    grp = gmsh.model.add_physical_group(2, surfaces)
    gmsh.model.set_physical_name(2, grp, "WALLS")


def main():
    """Run the example workflow."""
    try:
        # ---------------------------------------------------------------
        # Working directory (created next to this script)
        workdir = os.path.join(here, "abq")
        if not os.path.isdir(workdir):
            os.makedirs(workdir)
        os.chdir(workdir)

        # Material parameters (of the solid wall material)
        E = 210.0    # Young's modulus
        nu = 0.3     # Poisson's ratio

        # Honeycomb parameters
        h = 1.0             # vertical rib length
        l = 1.0             # inclined rib length
        theta_deg = -18.0   # re-entrant angle (< 0)
        depth = 1.0         # extrusion depth along y (the prism axis)
        thickness = 0.05    # shell wall thickness
        element_size = 0.25  # shell element size

        # ---------------------------------------------------------------
        # Build the model once (geometry, mesh, material, section, step). Only the
        # macroscopic strain changes between load cases, so the model is reused and
        # only the periodic boundary conditions are updated in the loop. The anchor
        # node used to suppress rigid-body motion is picked automatically.
        create_reentrant_honeycomb("RVE", h=h, l=l, theta_deg=theta_deg, depth=depth)
        unit_cell = GmshRVE("RVE")
        unit_cell.add_periodic_mesh_constraints()
        unit_cell.generate_mesh(element_size=element_size)

        unit_cell = AbaqusRVE(unit_cell)
        material = unit_cell.model.Material(name="LINEAR_ELASTICITY")
        material.Elastic(table=[[E, nu]])
        unit_cell.assign_shell_section("WALLS", "LINEAR_ELASTICITY", thickness)
        unit_cell.model.StaticStep(name="LOAD", previous="Initial")

        # Run the homogenization routine (six load cases in 3D)
        strain_magnitude = 0.03
        C = np.empty((6, 6))
        for i in range(6):
            print("\n------------------------------")
            print("Run load case {}.".format(i + 1))

            def strain_history(time, component=i):
                """Return the ramped macroscopic strain at ``time``."""
                strain = np.zeros(6)
                strain[component] = strain_magnitude
                return time*voigt_to_tensor(strain, "STRAIN")

            history = LoadHistory(np.array([0.0, 1.0]), strain_history)
            unit_cell.apply_periodic_BC("LOAD", history)

            unit_cell.run_analysis(field_outputs=["U"])
            stress = tensor_to_voigt(unit_cell.average_stress("LOAD"), "STRESS")
            C[:, i] = stress/strain_magnitude

            # Export the deformed walls to VTK
            if PARAQUS_INSTALLED:
                unit_cell.export_to_vtk(os.path.join(here, "vtk_output",
                                                     "case_{}".format(i + 1)),
                                        field_outputs=["U"], clear=True)


        if not PARAQUS_INSTALLED:
            print("Paraqus is not installed; VTK export was skipped.")

        # ---------------------------------------------------------------
        # Effective stiffness and engineering constants. The honeycomb lies in the
        # x-z plane and is extruded along y, so the response is orthotropic: auxetic
        # in the x-z plane and an ordinary (stiff) honeycomb prism along y.
        S = np.linalg.inv(C)
        Ex, Ez, Ey = 1.0/S[0, 0], 1.0/S[2, 2], 1.0/S[1, 1]
        nu_xz = -S[2, 0]/S[0, 0]     # z-strain per x-strain (in-plane)
        nu_zx = -S[0, 2]/S[2, 2]     # x-strain per z-strain (in-plane)

        print("\n------------------------------")
        print("Effective constitutive matrix (Voigt):")
        print(np.array2string(C, precision=4, suppress_small=True))
        print("  symmetry   max|C - C^T| = {:.2e}".format(np.max(np.abs(C - C.T))))
        print("\nIn-plane (x-z honeycomb) moduli and Poisson's ratios:")
        print("  E_x  = {:.4e}".format(Ex))
        print("  E_z  = {:.4e}".format(Ez))
        print("  nu_xz = {:+.4f}".format(nu_xz))
        print("  nu_zx = {:+.4f}".format(nu_zx))
        print("Out-of-plane (prism axis y):")
        print("  E_y  = {:.4e}   (stiff, non-auxetic)".format(Ey))
        if max(nu_xz, nu_zx) < 0.0:
            print("  --> negative in-plane Poisson's ratio: the honeycomb is AUXETIC.")

        # Export the effective stiffness
        result_path = os.path.join(here, "results")
        if not os.path.isdir(result_path):
            os.makedirs(result_path)
        export_to_csv_file(C, os.path.join(result_path, "auxetic_C.csv"))
    finally:
        if gmsh.is_initialized():
            gmsh.finalize()


if __name__ == "__main__":
    main()
