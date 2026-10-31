# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-only

"""Homogenize a continuum matrix reinforced by embedded beam fibers.

Two beam fibers span the cell in the z-direction and are coupled to the
continuum matrix by an embedded-region constraint. Their boundary nodes
inherit translational periodicity from the matrix, while their rotational
degrees of freedom remain part of the beam formulation. The example computes
the effective stiffness of the mixed-dimensional RVE.
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
from abqhom.utils import (export_to_csv_file, tensor_to_voigt,
                          voigt_to_tensor)


def create_geometry(model_name, edge_length=1.0, nodes_matrix=10,
                    nodes_fibers=3):
    """Create a hexahedral matrix containing two independent beam fibers."""

    # Initialize Gmsh and add a new model
    if not gmsh.is_initialized():
        gmsh.initialize()
    gmsh.model.add(model_name)

    gmsh.option.set_number("Geometry.Tolerance", 1e-5)

    x = y = z = 0.0  # Center point
    dx = dy = dz = edge_length  # Edge length

    # Add a cube to the model representing the matrix
    p1 = gmsh.model.geo.add_point(x - 0.5*dx, y - 0.5*dy, z - 0.5*dz)
    p2 = gmsh.model.geo.add_point(x - 0.5*dx, y - 0.5*dy, z + 0.5*dz)
    p3 = gmsh.model.geo.add_point(x + 0.5*dx, y - 0.5*dy, z + 0.5*dz)
    p4 = gmsh.model.geo.add_point(x + 0.5*dx, y - 0.5*dy, z - 0.5*dz)
    p5 = gmsh.model.geo.add_point(x - 0.5*dx, y + 0.5*dy, z - 0.5*dz)
    p6 = gmsh.model.geo.add_point(x - 0.5*dx, y + 0.5*dy, z + 0.5*dz)
    p7 = gmsh.model.geo.add_point(x + 0.5*dx, y + 0.5*dy, z + 0.5*dz)
    p8 = gmsh.model.geo.add_point(x + 0.5*dx, y + 0.5*dy, z - 0.5*dz)

    line_connectivity = [[p1, p2], [p2, p3], [p3, p4], [p4, p1],
                         [p5, p6], [p6, p7], [p7, p8], [p8, p5],
                         [p1, p5], [p2, p6], [p3, p7], [p4, p8]]
    lines = [gmsh.model.geo.add_line(*i) for i in line_connectivity]

    surfaces = []
    surface_connectivity = [[lines[0], lines[1], lines[2], lines[3]],
                            [lines[4], lines[5], lines[6], lines[7]],
                            [lines[0], lines[9], lines[4], lines[8]],
                            [lines[1], lines[10], lines[5], lines[9]],
                            [lines[2], lines[11], lines[6], lines[10]],
                            [lines[3], lines[8], lines[7], lines[11]]]
    for i in surface_connectivity:
        loop = gmsh.model.geo.add_curve_loop(i, reorient=True)
        surfaces.append(gmsh.model.geo.add_plane_surface([loop]))

    loop = gmsh.model.geo.add_surface_loop(surfaces)
    volume = gmsh.model.geo.add_volume([loop])

    # Add transfinite hex mesh constraints on the cube
    for i in lines:
        gmsh.model.geo.mesh.set_transfinite_curve(i, nodes_matrix)

    for i in surfaces:
        gmsh.model.geo.mesh.set_transfinite_surface(i)
        gmsh.model.geo.mesh.set_recombine(2, i)

    gmsh.model.geo.mesh.set_transfinite_volume(volume)

    # Add two fibers running the full height of the cell, so their end
    # nodes lie on the z-boundary faces (and, through the periodic tiling,
    # the fibers are continuous along z). The embedded-region constraint
    # eliminates those nodes' translational DOFs; apply_periodic_BC excludes
    # them from the translational equations, so the fibers may reach the
    # boundary.
    fiber_p1 = gmsh.model.geo.add_point(-0.25*dx, 0.25*dy, -0.5*dz)
    fiber_p2 = gmsh.model.geo.add_point(-0.25*dx, 0.25*dy, 0.5*dz)
    fiber_p3 = gmsh.model.geo.add_point(0.25*dx, -0.25*dy, -0.5*dz)
    fiber_p4 = gmsh.model.geo.add_point(0.25*dx, -0.25*dy, 0.5*dz)
    fiber_1 = gmsh.model.geo.add_line(fiber_p1, fiber_p2)
    fiber_2 = gmsh.model.geo.add_line(fiber_p3, fiber_p4)
    gmsh.model.geo.mesh.set_transfinite_curve(fiber_1, nodes_fibers)
    gmsh.model.geo.mesh.set_transfinite_curve(fiber_2, nodes_fibers)
    fibers = [fiber_1, fiber_2]

    # Synchronize the geo kernel to save the geometry
    gmsh.model.geo.synchronize()

    # Add physical groups
    matrix = gmsh.model.add_physical_group(3, [volume])
    gmsh.model.set_physical_name(3, matrix, "MATRIX")
    fibers = gmsh.model.add_physical_group(1, fibers)
    gmsh.model.set_physical_name(1, fibers, "FIBER")


def main():
    """Run the example workflow."""
    try:
        # ---------------------------------------------------------------
        # Working directory (created next to this script)
        workdir = os.path.join(here, "abq")
        if not os.path.isdir(workdir):
            os.makedirs(workdir)
        os.chdir(workdir)

        vtk_path = os.path.join(here, "vtk_output")

        # Material parameters: Young's modulus - Poisson's ratio
        parameters = {"MATRIX": [50.0, 0.35],
                      "FIBER":  [210, 0.2]}


        # ---------------------------------------------------------------
        # Build the model once (geometry, mesh, materials, sections, embedded
        # region, step). Only the macroscopic strain changes between load cases,
        # so the model is reused and only the periodic BCs are updated in the loop.
        create_geometry("RVE", nodes_matrix=9, nodes_fibers=5)
        unit_cell = GmshRVE("RVE")
        unit_cell.add_periodic_mesh_constraints()
        unit_cell.generate_mesh()

        unit_cell = AbaqusRVE(unit_cell)
        for name, (E, nu) in parameters.items():
            material = unit_cell.model.Material(name=name)
            material.Elastic(table=[[E, nu]])
        unit_cell.assign_solid_section("MATRIX", "MATRIX")
        unit_cell.assign_beam_section("FIBER", "FIBER", 0.1)

        # Embed the beam fibers into the continuum matrix. This constraint is kept
        # across load cases -- re-applying the periodic BC does not remove it.
        matrix_elements = unit_cell.instance.sets["ELSET_MATRIX"]
        fiber_elements = unit_cell.instance.sets["ELSET_FIBER"]
        unit_cell.model.EmbeddedRegion("COUPLING", fiber_elements, matrix_elements)
        unit_cell.model.StaticStep(name="LOAD", previous="Initial")

        # Run the six linear homogenization load cases.
        C = np.empty((6, 6))
        strain_magnitude = 0.03
        for i in range(6):
            print("\n------------------------------")
            print("Run case {}.".format(i + 1))

            def strain_history(time, component=i, magnitude=strain_magnitude):
                """Return the ramped macroscopic strain at ``time``."""
                strain = np.zeros(6)
                strain[component] = magnitude
                strain = voigt_to_tensor(strain, "STRAIN")
                return time*strain

            history = LoadHistory(np.array([0.0, 1.0]), strain_history)
            print("Strain = " + str(tensor_to_voigt(history(1.0), "STRAIN")
                                    .round(3).tolist()))

            unit_cell.apply_periodic_BC("LOAD", history)
            unit_cell.run_analysis(field_outputs=["U", "S"])

            stress = tensor_to_voigt(unit_cell.average_stress("LOAD"), "STRESS")
            C[:, i] = stress/strain_magnitude

            if PARAQUS_INSTALLED:
                unit_cell.export_to_vtk(
                    os.path.join(vtk_path, "case_{}".format(i + 1)),
                    field_outputs=["U", ("S", "MISES")],
                    element_sets=["MATRIX", "FIBER"], clear=True)


        if not PARAQUS_INSTALLED:
            print("Paraqus is not installed; VTK export was skipped.")

        print("\n------------------------------")
        print("Effective constitutive matrix (Voigt):")
        print(np.array2string(C, precision=4, suppress_small=True))

        # Store the effective stiffness
        result_dir = os.path.join(here, "results")
        if not os.path.isdir(result_dir):
            os.makedirs(result_dir)
        export_to_csv_file(C, os.path.join(result_dir, "stiffness.csv"))
    finally:
        if gmsh.is_initialized():
            gmsh.finalize()


if __name__ == "__main__":
    main()
