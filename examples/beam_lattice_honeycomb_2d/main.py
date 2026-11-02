# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-only

"""Homogenize a regular 2D beam honeycomb in the linear regime.

The struts intersect the cell boundary away from its corners. The example
demonstrates automatic boundary classification, periodic pairing, and
rigid-body anchoring for this general lattice topology. The resulting
effective stiffness is checked for in-plane isotropy.
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
from abqhom.gmshmodel import GmshRVE
from abqhom.utils import export_to_csv_file, tensor_to_voigt, voigt_to_tensor


def _clip_segment(p, q, box, tol=1e-9):
    """Clip a segment to a box with the Liang--Barsky method."""
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
    """Create a regular honeycomb beam lattice in a periodic cell.

    The infinite lattice is clipped to a rectangular period of dimensions
    ``3*a`` by ``sqrt(3)*a``. Opposite boundaries match without requiring a
    lattice node at a cell corner.
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
        """Return a coordinate key suitable for geometric deduplication."""
        return (round(pt[0], 7), round(pt[1], 7))

    point_tags = {}

    def get_point(pt):
        """Return the Gmsh tag of a deduplicated geometric point."""
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


def main():
    """Run the example workflow."""
    try:
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
        # Run the homogenization routine (three load cases in 2D: xx, yy, xy)
        components = [0, 1, 5]  # Voigt indices: xx, yy, xy(engineering)
        strain_magnitude = 0.03

        # Build the model once (geometry, mesh, material, section, step). Only the
        # macroscopic strain changes between load cases, so the model is reused and
        # only the periodic boundary conditions are updated in the loop. The anchor
        # is picked automatically -- no corner point required.
        create_honeycomb("RVE", a=a)
        unit_cell = GmshRVE("RVE")
        unit_cell.add_periodic_mesh_constraints()
        unit_cell.generate_mesh(element_size=element_size)

        unit_cell = AbaqusRVE(unit_cell)
        material = unit_cell.model.Material(name="LINEAR_ELASTICITY")
        material.Elastic(table=[[E, nu]])
        unit_cell.assign_beam_section("STRUTS", "LINEAR_ELASTICITY", radius)
        unit_cell.model.StaticStep(name="LOAD", previous="Initial")

        C = np.empty((3, 3))
        for col, comp in enumerate(components):
            print("\n------------------------------")
            print("Run load case {}.".format(col + 1))

            def strain_history(time, component=comp):
                """Return the ramped macroscopic strain at ``time``."""
                strain = np.zeros(6)
                strain[component] = strain_magnitude
                return time*voigt_to_tensor(strain, "STRAIN")

            history = LoadHistory(np.array([0.0, 1.0]), strain_history)
            unit_cell.apply_periodic_BC("LOAD", history)

            # Run the analysis and compute the average stress
            unit_cell.run_analysis(field_outputs=["U"])
            stress = unit_cell.average_stress("LOAD")
            stress = tensor_to_voigt(stress, "STRESS")
            print("Stress = " + str(stress.round(8).tolist()))

            C[:, col] = stress/strain_magnitude

            # Export the deformed struts to VTK (the strut radius is taken
            # automatically from the assigned beam section)
            unit_cell.export_to_vtk(os.path.join(here, "vtk_output",
                                                 "case_{}".format(col + 1)),
                                    field_outputs=["U"], clear=True)

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
    finally:
        if gmsh.is_initialized():
            gmsh.finalize()


if __name__ == "__main__":
    main()
