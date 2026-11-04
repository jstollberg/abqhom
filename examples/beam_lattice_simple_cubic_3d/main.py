# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-or-later

"""Homogenize a simple-cubic beam lattice in the linear regime.

The twelve edge struts are shared by four neighboring cells, making this the
simplest example of the automatic boundary-element reduction. The numerical
stiffness is compared with ``C11 = E*pi*r**2/L**2``, while the transverse
coupling terms are expected to vanish.
"""

import os
import sys

import gmsh
import numpy as np

from abqhom.abqmodel import AbaqusRVE, LoadHistory
from abqhom.gmshmodel import GmshRVE
from abqhom.utils import export_to_csv_file, tensor_to_voigt, voigt_to_tensor

here = os.getcwd()


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
    corners = [g.add_point(x, y, z)
               for x in (-h, h) for y in (-h, h) for z in (-h, h)]

    # Connect the corners that differ in exactly one coordinate (edges)
    struts = []
    coords = [(x, y, z) for x in (-h, h) for y in (-h, h) for z in (-h, h)]
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
    gmsh.model.set_physical_name(1, group, "STRUTS")


def main():
    """Run the example workflow."""
    try:
        # ---------------------------------------------------------------
        # Working directory (created next to this script)
        workdir = os.path.join(here, "abq")
        if not os.path.isdir(workdir):
            os.makedirs(workdir)
        os.chdir(workdir)

        # Material and geometry
        E = 210.0        # Young's modulus of the strut material
        nu = 0.3         # Poisson's ratio
        length = 10.0    # edge length of the cubic unit cell
        radius = 0.5     # strut radius
        element_size = 1.0

        # ---------------------------------------------------------------
        # Build the model once and update only the periodic boundary conditions
        # between load cases.
        model_name = "RVE"
        create_geometry(model_name, length)
        unit_cell = GmshRVE(model_name)
        unit_cell.add_periodic_mesh_constraints()
        unit_cell.generate_mesh(element_size=element_size)

        unit_cell = AbaqusRVE(unit_cell)
        material = unit_cell.model.Material(name="STRUT")
        material.Elastic(table=[[E, nu]])
        unit_cell.assign_beam_section("STRUTS", "STRUT", radius)
        unit_cell.model.StaticStep(name="LOAD", previous="Initial")

        # Run the homogenization routine (six load cases in 3D)
        strain_magnitude = 0.03
        component_names = ("xx", "yy", "zz", "yz", "xz", "xy")
        C = np.zeros((6, 6))
        for i, component_name in enumerate(component_names):
            print("\n------------------------------")
            print("Run strain case {} ({}).".format(i + 1, component_name))

            def strain_history(time, component=i):
                """Return the ramped macroscopic strain at ``time``."""
                strain = np.zeros(6)
                strain[component] = strain_magnitude
                return time*voigt_to_tensor(strain, "STRAIN")

            history = LoadHistory(np.array([0.0, 1.0]), strain_history)
            unit_cell.apply_periodic_bc("LOAD", history)
            print("  Run Abaqus analysis.")
            sys.stdout.flush()
            unit_cell.run_analysis(field_outputs=["U"])

            stress = tensor_to_voigt(
                unit_cell.average_stress("LOAD"), "STRESS")
            C[:, i] = stress/strain_magnitude

            print("  Export VTK results.")
            sys.stdout.flush()
            unit_cell.export_to_vtk(os.path.join(here, "vtk_output",
                                                 component_name),
                                    field_outputs=["U"], clear=True)

        # ---------------------------------------------------------------
        C11_analytical = E*np.pi*radius**2/length**2
        print("\n------------------------------")
        print("Effective constitutive matrix (Voigt):")
        print(np.array2string(C, precision=4, suppress_small=True))
        print("\nComparison with the analytical simple-cubic result:")
        print("  C11 (FE)         = {:.6e}".format(C[0, 0]))
        print("  C11 (analytical) = {:.6e}".format(C11_analytical))
        print("  relative error   = {:.2e}".format(
            abs(C[0, 0] - C11_analytical)/C11_analytical))
        print("  max |C12,C13,C23| (should vanish) = {:.2e}".format(
            max(abs(C[0, 1]), abs(C[0, 2]), abs(C[1, 2]))))

        # Store the effective stiffness
        result_dir = os.path.join(here, "results")
        if not os.path.isdir(result_dir):
            os.makedirs(result_dir)
        export_to_csv_file(C, os.path.join(result_dir, "stiffness.csv"))
    finally:
        if gmsh.is_initialized():
            gmsh.finalize()


if __name__ == "__main__":
    # Abaqus replaces stdout; restore the shell that launched the example.
    if sys.__stdout__ is not None:
        sys.stdout = sys.__stdout__
    main()
