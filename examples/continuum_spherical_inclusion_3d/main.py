# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-or-later

"""Homogenize a two-phase spherical-inclusion cell at finite strain.

Displacement and stress control are applied successively to verify that the
same effective deformation and stress state is recovered.
"""

import os
import sys

import gmsh
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from abaqusConstants import ISOTROPIC, NEO_HOOKE, OFF, ON, VOLUMETRIC_DATA

from abqhom.abqmodel import AbaqusRVE, LoadHistory
from abqhom.gmshmodel import GmshRVE
from abqhom.utils import export_to_csv_file

here = os.getcwd()


def create_geometry(model_name, dx=1.0, dy=1.0, dz=1.0,
                    inclusion_fraction=0.2):
    """Create a cuboidal matrix with a centered spherical inclusion."""
    # Determine the sphere radius from the prescribed inclusion volume
    # fraction. The cell dimensions may differ, although this example uses a
    # cube below.
    radius = ((0.75*inclusion_fraction*dx*dy*dz)/np.pi)**(1.0/3.0)
    if 2*radius >= min(dx, dy, dz):
        raise ValueError("Inclusion is bigger than the RVE.")

    if not gmsh.is_initialized():
        gmsh.initialize()
    if model_name in gmsh.model.list():
        gmsh.model.set_current(model_name)
        gmsh.model.remove()
    gmsh.model.add(model_name)

    gmsh.option.set_number("Geometry.Tolerance", 1e-5)

    x = y = z = 0.0

    box = gmsh.model.occ.add_box(
        x - dx/2.0, y - dy/2.0, z - dz/2.0, dx, dy, dz)
    inclusion = gmsh.model.occ.add_sphere(x, y, z, radius)

    # Subtract a copy of the sphere from the box while retaining the original
    # sphere. The two resulting volumes share a conforming interface.
    matrix, _ = gmsh.model.occ.cut(
        [(3, box)], [(3, inclusion)], removeTool=False)
    gmsh.model.occ.synchronize()

    # Physical names become Abaqus element sets and are used below to assign
    # the two constituent materials.
    group = gmsh.model.add_physical_group(3, [tag for _, tag in matrix])
    gmsh.model.set_physical_name(3, group, "MATRIX")
    group = gmsh.model.add_physical_group(3, [inclusion])
    gmsh.model.set_physical_name(3, group, "INCLUSION")


def evaluate_response(unit_cell, step_name):
    """Return the time, deformation-gradient, and stress histories."""

    # Keep the ODB open while both effective quantities are evaluated at
    # every frame. The post-processing methods then reuse the same handle.
    # At finite strain, average_stress returns the effective first
    # Piola--Kirchhoff stress that is work-conjugate to F.
    unit_cell.open_odb_file()
    try:
        frames = unit_cell.odb.steps[step_name].frames
        times = []
        deformation_gradients = []
        stresses = []
        for index, frame in enumerate(frames):
            times.append(frame.frameValue)
            deformation_gradients.append(
                unit_cell.average_deformation_gradient(
                    step_name, frame=index))
            stresses.append(unit_cell.average_stress(step_name, frame=index))
    finally:
        unit_cell.close_odb_file()
    return (np.asarray(times), np.asarray(deformation_gradients),
            np.asarray(stresses))


def main():
    """Run the example workflow."""
    try:
        # ---------------------------------------------------------------
        # 1. Define the analysis and create its working directory.
        # ---------------------------------------------------------------
        workdir = os.path.join(here, "abq")
        if not os.path.isdir(workdir):
            os.makedirs(workdir)
        os.chdir(workdir)

        # RVE and material parameters (E, nu).
        dx = dy = dz = 10.0
        inclusion_fraction = 0.3       # volume fraction of the inclusion
        parameters = {"MATRIX":    [3.5, 0.35],    # soft matrix
                      "INCLUSION": [230.0, 0.2]}   # stiff inclusion
        element_size = 1.0
        final_displacement_gradient = 0.2

        # ---------------------------------------------------------------
        # 2. Build the geometry and impose a matching periodic mesh.
        # ---------------------------------------------------------------
        # Matching nodes on opposing faces permit direct periodic equation
        # constraints after the Gmsh mesh has been imported into Abaqus.
        create_geometry("RVE", dx=dx, dy=dy, dz=dz,
                        inclusion_fraction=inclusion_fraction)
        unit_cell = GmshRVE("RVE")
        unit_cell.add_periodic_mesh_constraints()
        unit_cell.generate_mesh(element_size=element_size)

        # ---------------------------------------------------------------
        # 3. Create the Abaqus model, materials, and nonlinear load step.
        # ---------------------------------------------------------------
        # Transfer the mesh and assign a compressible Neo-Hookean material to
        # each constituent. Abaqus uses C10 = G/2 and D1 = 2/K for the chosen
        # material form and performs the constitutive integration.
        unit_cell = AbaqusRVE(unit_cell)
        for name, (E, nu) in parameters.items():
            shear_modulus = E/(2.0*(1.0 + nu))
            bulk_modulus = E/(3.0*(1.0 - 2.0*nu))
            material = unit_cell.model.Material(name=name)
            material.Hyperelastic(
                materialType=ISOTROPIC,
                table=((0.5*shear_modulus, 2.0/bulk_modulus),),
                testData=OFF, type=NEO_HOOKE,
                volumetricResponse=VOLUMETRIC_DATA)
            unit_cell.assign_solid_section(name, name)

        step_name = "LOAD"
        increment = 0.01

        # Equal minimum, initial, and maximum increments give the two control
        # runs the same pseudo-time discretization for a framewise comparison.
        unit_cell.model.StaticStep(
            name=step_name, previous="Initial", nlgeom=ON, maxNumInc=200,
            initialInc=increment, minInc=increment, maxInc=increment)

        # ---------------------------------------------------------------
        # 4. Generate a reference solution under deformation control.
        # ---------------------------------------------------------------
        # The displacement-controlled API receives the macroscopic
        # displacement gradient H = F - I. Only H11 is increased here, while
        # all other components remain fixed at zero. The complete effective
        # stress tensor is retained at every converged frame.
        displacement_history = LoadHistory(
            np.array([0.0, 1.0]),
            lambda pseudo_time: np.diag(
                [final_displacement_gradient*pseudo_time, 0.0, 0.0]))
        unit_cell.apply_periodic_displacement_bc(
            step_name, displacement_history)
        print("\nRun displacement-controlled finite-strain load path.")
        print("  Run Abaqus analysis.")
        sys.stdout.flush()
        unit_cell.run_analysis(field_outputs=["U", "S"])

        print("  Post-process displacement-controlled response.")
        sys.stdout.flush()
        displacement_times, target_deformations, target_stresses = (
            evaluate_response(unit_cell, step_name))

        # ---------------------------------------------------------------
        # 5. Recover the same path under first Piola--Kirchhoff stress control.
        # ---------------------------------------------------------------
        # Use the complete nonlinear stress history of the reference path.
        # Scaling only its final value would generally define a different
        # multiaxial path and would not reproduce the intermediate states.
        def stress_path(pseudo_time):
            """Interpolate the displacement-controlled effective stress."""

            # Abaqus amplitudes are scalar, so interpolate every tensor
            # component independently at the requested pseudo-time.
            values = [
                np.interp(pseudo_time, displacement_times,
                          target_stresses[:, row, column])
                for row in range(3) for column in range(3)
            ]
            return np.asarray(values).reshape((3, 3))

        stress_history = LoadHistory(
            displacement_times, stress_path)

        # The periodic traction conditions replace the preceding macroscopic
        # displacement constraints. They create generalized loads conjugate
        # to the independent components of F, while the microscopic periodic
        # equations retain their kinematic form.
        unit_cell.apply_periodic_traction_bc(step_name, stress_history)
        print("\nRun stress-controlled finite-strain load path.")
        print("  Run Abaqus analysis.")
        sys.stdout.flush()
        unit_cell.run_analysis(field_outputs=["U", "S"])

        print("  Post-process stress-controlled response.")
        sys.stdout.flush()
        _, recovered_deformations, recovered_stresses = (
            evaluate_response(unit_cell, step_name))

        # ---------------------------------------------------------------
        # 6. Compare, report, and export both solutions.
        # ---------------------------------------------------------------
        # Compare the complete tensor histories rather than only their final
        # states. One global scale per history avoids division by the zero
        # initial displacement gradient and stress.
        displacement_curve = np.column_stack(
            (target_deformations[:, 0, 0], target_stresses[:, 0, 0]))
        stress_curve = np.column_stack(
            (recovered_deformations[:, 0, 0],
             recovered_stresses[:, 0, 0]))

        deformation_error = max(
            np.linalg.norm(recovered - target)
            for recovered, target in zip(
                recovered_deformations, target_deformations))
        deformation_scale = max(
            np.linalg.norm(deformation - np.eye(3))
            for deformation in target_deformations)
        deformation_error /= deformation_scale

        stress_error = max(
            np.linalg.norm(recovered - target)
            for recovered, target in zip(recovered_stresses, target_stresses))
        stress_scale = max(
            np.linalg.norm(stress) for stress in target_stresses)
        stress_error /= stress_scale

        # The plotted axial components provide a compact visual comparison.
        # The residuals below nevertheless include all tensor components.
        print("\n------------------------------")
        print("Finite-strain stress-control verification:")
        print("  displacement-control frames = {}".format(
            len(displacement_curve)))
        print("  stress-control frames       = {}".format(len(stress_curve)))
        print("  target F11                  = {:.6f}".format(
            target_deformations[-1, 0, 0]))
        print("  recovered F11               = {:.6f}".format(
            recovered_deformations[-1, 0, 0]))
        print("  target P11                  = {:.6e}".format(
            target_stresses[-1, 0, 0]))
        print("  recovered P11               = {:.6e}".format(
            recovered_stresses[-1, 0, 0]))
        print("  maximum relative deformation error = {:.2e}".format(
            deformation_error))
        print("  maximum relative stress error      = {:.2e}".format(
            stress_error))

        result_dir = os.path.join(here, "results")
        if not os.path.isdir(result_dir):
            os.makedirs(result_dir)
        export_to_csv_file(
            displacement_curve,
            os.path.join(result_dir, "response_displacement.csv"),
            header="F_11,P_11")
        export_to_csv_file(
            stress_curve, os.path.join(result_dir, "response_stress.csv"),
            header="F_11,P_11")
        export_to_csv_file(
            np.array([[deformation_error, stress_error]]),
            os.path.join(result_dir, "verification.csv"),
            header=("maximum_relative_deformation_error,"
                    "maximum_relative_stress_error"))

        # Overlay the reference and recovered axial response curves.
        fig, ax = plt.subplots(figsize=(6.0, 5.0))
        ax.plot(displacement_curve[:, 0], displacement_curve[:, 1], "-",
                color="#1f6feb", lw=1.4, label="Displacement control")
        ax.plot(stress_curve[:, 0], stress_curve[:, 1], "--",
                color="#d97706", lw=1.4, label="Stress control")
        ax.axhline(0.0, color="0.6", lw=0.6)
        ax.axvline(1.0, color="0.6", lw=0.6)
        ax.set_xlabel(r"Macroscopic deformation gradient  $\bar{F}_{11}$")
        ax.set_ylabel(
            r"Effective first Piola-Kirchhoff stress  $\bar{P}_{11}$")
        ax.set_title("Finite-strain response of a spherical-inclusion cell")
        ax.grid(True, alpha=0.3)
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(os.path.join(result_dir, "response.png"), dpi=150)
        plt.close(fig)

        # Export the final stress-controlled analysis for visualization.
        print("  Export VTK results.")
        sys.stdout.flush()
        unit_cell.export_to_vtk(
            os.path.join(here, "vtk_output"),
            field_outputs=["U", ("S", "MISES")],
            element_sets=["MATRIX", "INCLUSION"], clear=True)

        print("*** Finished ***")
    finally:
        if gmsh.is_initialized():
            gmsh.finalize()


if __name__ == "__main__":
    # Abaqus replaces stdout; restore the shell that launched the example.
    if sys.__stdout__ is not None:
        sys.stdout = sys.__stdout__
    main()
