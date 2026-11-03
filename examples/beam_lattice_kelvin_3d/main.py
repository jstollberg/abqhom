# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-or-later

"""Homogenize a Kelvin beam lattice at finite strain.

The vertices of the truncated octahedron lie on the faces of the cuboidal
cell rather than at its corners. The example demonstrates automatic boundary
classification and periodic pairing for this topology and computes a
monotonic finite-strain response.
"""

import os
from itertools import permutations, product

import gmsh
import numpy as np
from abaqusConstants import ISOTROPIC, NEO_HOOKE, OFF, ON, VOLUMETRIC_DATA

from abqhom.abqmodel import AbaqusRVE, LoadHistory
from abqhom.gmshmodel import GmshRVE
from abqhom.utils import export_to_csv_file

here = os.getcwd()


def create_kelvin(model_name, edge=4.0):
    """Create a Kelvin beam lattice in a cuboidal periodic cell.

    The 24 vertices are the signed permutations of
    ``(0, 1, 2)*edge/4``. Struts connect vertices separated by
    ``sqrt(2)*edge/4``.
    """
    if not gmsh.is_initialized():
        gmsh.initialize()
    if model_name in gmsh.model.list():
        gmsh.model.set_current(model_name)
        gmsh.model.remove()
    gmsh.model.add(model_name)
    g = gmsh.model.geo

    s = edge/4.0
    verts = set()
    for perm in set(permutations((0, 1, 2))):
        for signs in product((1, -1), repeat=3):
            verts.add(tuple(round(perm[i]*signs[i]*s, 9) for i in range(3)))
    verts = sorted(verts)

    tag = {v: g.add_point(v[0], v[1], v[2]) for v in verts}
    target = np.sqrt(2.0)*s
    lines = []
    for i in range(len(verts)):
        for j in range(i + 1, len(verts)):
            d = np.linalg.norm(np.array(verts[i]) - np.array(verts[j]))
            if abs(d - target) < 1e-6:
                lines.append(g.add_line(tag[verts[i]], tag[verts[j]]))

    g.synchronize()
    grp = gmsh.model.add_physical_group(1, lines)
    gmsh.model.set_physical_name(1, grp, "STRUTS")


def get_frame_times(unit_cell, step_name):
    """Return the output-frame times and close the ODB."""
    unit_cell.open_odb_file()
    try:
        return [frame.frameValue
                for frame in unit_cell.odb.steps[step_name].frames]
    finally:
        unit_cell.close_odb_file()


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
        final_displacement_gradient = 0.2

        # Lattice parameters
        edge = 4.0          # cubic unit-cell edge length
        radius = 0.05       # strut cross-section radius
        element_size = 0.6  # beam element length

        # ---------------------------------------------------------------
        # The anchor is selected automatically, so this face-vertex lattice
        # requires no corner node.
        create_kelvin("RVE", edge=edge)
        unit_cell = GmshRVE("RVE")
        unit_cell.add_periodic_mesh_constraints()
        unit_cell.generate_mesh(element_size=element_size)

        unit_cell = AbaqusRVE(unit_cell, custom_types={1: "B31"})
        shear_modulus = E/(2.0*(1.0 + nu))
        bulk_modulus = E/(3.0*(1.0 - 2.0*nu))
        material = unit_cell.model.Material(name="NEO_HOOKE")
        material.Hyperelastic(
            materialType=ISOTROPIC,
            table=((0.5*shear_modulus, 2.0/bulk_modulus),),
            testData=OFF, type=NEO_HOOKE,
            volumetricResponse=VOLUMETRIC_DATA)
        unit_cell.assign_beam_section("STRUTS", "NEO_HOOKE", radius)

        step_name = "LOAD"
        increment = 0.01
        unit_cell.model.StaticStep(
            name=step_name, previous="Initial", nlgeom=ON, maxNumInc=200,
            initialInc=increment, minInc=increment, maxInc=increment)

        history = LoadHistory(
            np.array([0.0, 1.0]),
            lambda pseudo_time: np.diag(
                [final_displacement_gradient*pseudo_time, 0.0, 0.0]))
        unit_cell.apply_periodic_bc(step_name, history)
        unit_cell.run_analysis(field_outputs=["U", "S"])

        frame_times = get_frame_times(unit_cell, step_name)

        curve = []
        for index, frame_time in enumerate(frame_times):
            F11 = 1.0 + final_displacement_gradient*frame_time
            stress = unit_cell.average_stress(step_name, frame=index)
            curve.append([F11, stress[0, 0]])
        curve = np.asarray(curve)

        print("\n------------------------------")
        print("Finite-strain response ({} frames):".format(len(curve)))
        print("  final F11 = {:.4f}".format(curve[-1, 0]))
        print("  final P11 = {:.6e}".format(curve[-1, 1]))

        result_dir = os.path.join(here, "results")
        if not os.path.isdir(result_dir):
            os.makedirs(result_dir)
        export_to_csv_file(curve, os.path.join(result_dir, "response.csv"))

        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(6.0, 5.0))
        ax.plot(curve[:, 0], curve[:, 1], "-", color="#1f6feb", lw=1.4)
        ax.axhline(0.0, color="0.6", lw=0.6)
        ax.axvline(1.0, color="0.6", lw=0.6)
        ax.set_xlabel(r"Macroscopic deformation gradient  $\bar{F}_{11}$")
        ax.set_ylabel(
            r"Effective first Piola-Kirchhoff stress  $\bar{P}_{11}$")
        ax.set_title("Finite-strain response of a Kelvin lattice")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(os.path.join(result_dir, "response.png"), dpi=150)
        plt.close(fig)

        unit_cell.export_to_vtk(
            os.path.join(here, "vtk_output"),
            field_outputs=["U", ("S", "MISES")], clear=True)

        print("*** Finished ***")
    finally:
        if gmsh.is_initialized():
            gmsh.finalize()


if __name__ == "__main__":
    main()
