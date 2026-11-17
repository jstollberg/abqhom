# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-or-later

"""Homogenize a complete simple-cubic beam lattice in the linear regime.

The twelve edge struts are shared by four neighboring cells, making this the
simplest example of ownership weighting. The complete cell is retained and
the Young's modulus of every strut is scaled by its ownership factor 1/4.
The numerical stiffness is compared with ``C11 = E*pi*r**2/L**2``, while the
transverse coupling terms are expected to vanish. A subsequent uniaxial
stress-controlled case verifies the recovered effective stress and strain.
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
        edge_ownership = 0.25  # four neighboring cells share each edge

        # ---------------------------------------------------------------
        # Build the model once and update only the periodic boundary conditions
        # between load cases.
        model_name = "RVE"
        create_geometry(model_name, length)
        unit_cell = GmshRVE(model_name)
        unit_cell.add_periodic_mesh_constraints()
        unit_cell.generate_mesh(element_size=element_size)

        # Retain the complete periodic representation. Each edge strut is
        # shared by four cells, so scaling E by 1/4 weights its complete
        # linear-elastic beam response by the corresponding ownership factor.
        unit_cell = AbaqusRVE(unit_cell, remove_secondary=False)
        material = unit_cell.model.Material(name="WEIGHTED_STRUT")
        material.Elastic(table=[[edge_ownership*E, nu]])
        unit_cell.assign_beam_section(
            "STRUTS", "WEIGHTED_STRUT", radius)
        unit_cell.model.StaticStep(name="LOAD", previous="Initial")

        # Run the six independent three-dimensional strain cases. Each
        # recovered stress vector forms one column of the effective stiffness
        # matrix in Voigt notation.
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
            unit_cell.apply_periodic_displacement_bc("LOAD", history)
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

        # Verify the inverse response by prescribing the analytical axial
        # stress corresponding to the strain used above. The remaining free
        # macroscopic controls impose zero transverse and shear stress, and
        # the solution should recover the original axial strain.
        stress_magnitude = C11_analytical*strain_magnitude
        print("\n------------------------------")
        print("Run stress-controlled case (xx).")

        def stress_history(time):
            """Return the ramped macroscopic Cauchy stress."""
            stress = np.zeros((3, 3))
            stress[0, 0] = time*stress_magnitude
            return stress

        history = LoadHistory(np.array([0.0, 1.0]), stress_history)
        unit_cell.apply_periodic_stress_bc("LOAD", history)
        print("  Run Abaqus analysis.")
        sys.stdout.flush()
        unit_cell.run_analysis(field_outputs=["U"])

        effective_stress = unit_cell.average_stress("LOAD")
        deformation_gradient = unit_cell.average_deformation_gradient(
            "LOAD")
        recovered_strain = deformation_gradient[0, 0] - 1.0
        stress_residual = (
            np.linalg.norm(effective_stress - stress_history(1.0))
            / stress_magnitude)
        strain_error = (
            abs(recovered_strain - strain_magnitude)/strain_magnitude)

        print("  Export VTK results.")
        sys.stdout.flush()
        unit_cell.export_to_vtk(
            os.path.join(here, "vtk_output", "stress_xx"),
            field_outputs=["U"], clear=True)

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
        print("\nPeriodic stress-control verification:")
        print("  prescribed sigma11 = {:.6e}".format(stress_magnitude))
        print("  recovered sigma11  = {:.6e}".format(
            effective_stress[0, 0]))
        print("  recovered epsilon11 = {:.6e}".format(recovered_strain))
        print("  relative stress residual = {:.2e}".format(stress_residual))
        print("  relative strain error    = {:.2e}".format(strain_error))

        # Store the effective stiffness and stress-control verification.
        result_dir = os.path.join(here, "results")
        if not os.path.isdir(result_dir):
            os.makedirs(result_dir)
        export_to_csv_file(C, os.path.join(result_dir, "stiffness.csv"))
        stress_control = np.array([[stress_magnitude,
                                    effective_stress[0, 0],
                                    strain_magnitude,
                                    recovered_strain,
                                    stress_residual,
                                    strain_error]])
        export_to_csv_file(
            stress_control, os.path.join(result_dir, "stress_control.csv"),
            header=("prescribed_sigma_11,recovered_sigma_11,"
                    "expected_epsilon_11,recovered_epsilon_11,"
                    "stress_residual,strain_error"))
    finally:
        if gmsh.is_initialized():
            gmsh.finalize()


if __name__ == "__main__":
    # Abaqus replaces stdout; restore the shell that launched the example.
    if sys.__stdout__ is not None:
        sys.stdout = sys.__stdout__
    main()
