# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-only

"""Homogenize a two-phase plane-strain composite in the linear regime.

The cell contains a square arrangement of stiff circular inclusions in a soft
matrix. The example computes the effective in-plane stiffness and verifies
the Voigt and Reuss matrix bounds.
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
from abqhom.utils import (export_to_csv_file, reuss, tensor_to_voigt, voigt,
                          voigt_to_tensor)


def create_geometry(model_name, dx, dy, fiber_fraction):
    """Create a square array of circular inclusions."""
    radius = np.sqrt((fiber_fraction*dx*dy)/(2*np.pi))
    if radius >= 0.5*min(dx, dy):
        raise ValueError("Fibers overlap.")

    if not gmsh.is_initialized():
        gmsh.initialize()
    gmsh.model.add(model_name)
    gmsh.option.set_number("Geometry.Tolerance", 1e-5)

    rectangle = gmsh.model.occ.add_rectangle(-0.5*dx, -0.5*dy, 0.0, dx, dy)
    disks = [gmsh.model.occ.add_disk(cx, cy, 0.0, radius, radius)
             for cx, cy in [(0.0, 0.0),
                            (-0.5*dx, -0.5*dy), (0.5*dx, -0.5*dy),
                            (0.5*dx, 0.5*dy), (-0.5*dx, 0.5*dy)]]

    fibers, _ = gmsh.model.occ.intersect(
        [(2, d) for d in disks], [(2, rectangle)], removeTool=False)
    matrix, _ = gmsh.model.occ.cut([(2, rectangle)], fibers, removeTool=False)
    gmsh.model.occ.synchronize()

    group = gmsh.model.add_physical_group(2, [i for _, i in matrix])
    gmsh.model.set_physical_name(2, group, "MATRIX")
    group = gmsh.model.add_physical_group(2, [i for _, i in fibers])
    gmsh.model.set_physical_name(2, group, "FIBER")


def isotropic_plane_strain_stiffness(E, nu):
    """Return the isotropic plane-strain stiffness in Voigt notation."""
    C = np.array([[1 - nu, nu, 0.0],
                  [nu, 1 - nu, 0.0],
                  [0.0, 0.0, 0.5*(1 - 2*nu)]])
    return C*E/((1 + nu)*(1 - 2*nu))


def main():
    """Run the example workflow."""
    try:
        # ---------------------------------------------------------------
        # Working directory (created next to this script)
        workdir = os.path.join(here, "abq")
        if not os.path.isdir(workdir):
            os.makedirs(workdir)
        os.chdir(workdir)

        # RVE and material parameters of the isotropic phases (E, nu)
        dx = dy = 10.0
        fiber_fraction = 0.4          # change this (or loop) for other fractions
        parameters = {"MATRIX": [3.5, 0.35],    # soft matrix
                      "FIBER":  [230.0, 0.2]}   # stiff inclusions

        # ---------------------------------------------------------------
        # Run the homogenization routine (three load cases in 2D: xx, yy, xy)
        components = [0, 1, 5]  # Voigt indices: xx, yy, xy(engineering)
        strain_magnitude = 0.03

        # Build the model once (geometry, mesh, materials, sections, step). Only the
        # macroscopic strain changes between load cases, so the model is reused and
        # only the periodic boundary conditions are updated in the loop.
        create_geometry("RVE", dx, dy, fiber_fraction)
        unit_cell = GmshRVE("RVE")
        unit_cell.add_periodic_mesh_constraints()
        unit_cell.generate_mesh(element_size=0.4)

        unit_cell = AbaqusRVE(unit_cell)
        for name, (E, nu) in parameters.items():
            material = unit_cell.model.Material(name=name)
            material.Elastic(table=[[E, nu]])
            unit_cell.assign_solid_section(name, name)
        unit_cell.model.StaticStep(name="LOAD", previous="Initial")

        C = np.zeros((3, 3))
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
            unit_cell.run_analysis(field_outputs=["U", "S"])

            stress = tensor_to_voigt(unit_cell.average_stress("LOAD"), "STRESS")
            C[:, col] = stress/strain_magnitude

            if PARAQUS_INSTALLED:
                unit_cell.export_to_vtk(
                    os.path.join(here, "vtk_output", "case_{}".format(col + 1)),
                    field_outputs=["U", ("S", "MISES")],
                    element_sets=["MATRIX", "FIBER"], clear=True)


        if not PARAQUS_INSTALLED:
            print("Paraqus is not installed; VTK export was skipped.")

        # ---------------------------------------------------------------
        # Voigt / Reuss bounds of the two-phase material
        C_matrix = isotropic_plane_strain_stiffness(*parameters["MATRIX"])
        C_fiber = isotropic_plane_strain_stiffness(*parameters["FIBER"])
        C_reuss = reuss(fiber_fraction, C_fiber, C_matrix)
        C_voigt = voigt(fiber_fraction, C_fiber, C_matrix)

        C_symmetric = 0.5*(C + C.T)
        lower_margin = np.min(np.linalg.eigvalsh(C_symmetric - C_reuss))
        upper_margin = np.min(np.linalg.eigvalsh(C_voigt - C_symmetric))
        print("\n------------------------------")
        print("Effective constitutive matrix (Voigt):")
        print(np.array2string(C, precision=4, suppress_small=True))
        print("  symmetry   max|C - C^T| = {:.2e}".format(
            np.max(np.abs(C - C.T))))
        print("\nVoigt--Reuss matrix-bound margins (must be non-negative):")
        print("  min eig(C - C_Reuss) = {:.6e}".format(lower_margin))
        print("  min eig(C_Voigt - C) = {:.6e}".format(upper_margin))

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
