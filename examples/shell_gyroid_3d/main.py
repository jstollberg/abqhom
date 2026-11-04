# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-or-later

"""Homogenize an imported sheet-gyroid shell at finite strain.

The periodic surface mesh is generated separately by ``generate_mesh.py``.
This example imports the mesh through ``GmshRVE.read_from_mesh``, verifies its
periodicity, and computes the response along a monotonic loading path.
"""

import os
import sys

import gmsh
import numpy as np
from abaqusConstants import ISOTROPIC, NEO_HOOKE, OFF, ON, VOLUMETRIC_DATA

from abqhom.abqmodel import AbaqusRVE, LoadHistory
from abqhom.gmshmodel import GmshRVE
from abqhom.utils import export_to_csv_file

here = os.getcwd()


def main():
    """Run the example workflow."""
    try:
        # ---------------------------------------------------------------
        # Working directory (created next to this script)
        workdir = os.path.join(here, "abq")
        if not os.path.isdir(workdir):
            os.makedirs(workdir)
        os.chdir(workdir)

        # Shell thickness and material
        thickness = 0.05     # shell thickness
        E = 210.0            # Young's modulus of the shell material
        nu = 0.3             # Poisson's ratio
        final_displacement_gradient = 0.2

        # Load the pre-computed microgen mesh (physical group 'GYROID').
        mesh_file = os.path.join(here, "gyroid.msh")
        if not os.path.isfile(mesh_file):
            raise RuntimeError(
                "The mesh file '{}' does not exist. Generate it first by "
                "running generate_mesh.py in an environment that has microgen "
                "installed."
                .format(mesh_file))

        gmsh.initialize()

        # Merge the mesh into the common model name used by all examples.
        model_name = "RVE"
        gmsh.model.add(model_name)
        gmsh.merge(mesh_file)
        unit_cell = GmshRVE.read_from_mesh(model_name)
        if not unit_cell.periodic_mesh:
            raise RuntimeError("The imported mesh is not periodic.")

        # ---------------------------------------------------------------
        # Build the Abaqus model and assign the shell material and section.
        cell = AbaqusRVE(unit_cell)
        shear_modulus = E/(2.0*(1.0 + nu))
        bulk_modulus = E/(3.0*(1.0 - 2.0*nu))
        material = cell.model.Material(name="SHELL")
        material.Hyperelastic(
            materialType=ISOTROPIC,
            table=((0.5*shear_modulus, 2.0/bulk_modulus),),
            testData=OFF, type=NEO_HOOKE,
            volumetricResponse=VOLUMETRIC_DATA)
        cell.assign_shell_section("GYROID", "SHELL", thickness)

        step_name = "LOAD"
        increment = 0.01
        cell.model.StaticStep(
            name=step_name, previous="Initial", nlgeom=ON, maxNumInc=200,
            initialInc=increment, minInc=increment, maxInc=increment)

        history = LoadHistory(
            np.array([0.0, 1.0]),
            lambda pseudo_time: np.diag(
                [final_displacement_gradient*pseudo_time, 0.0, 0.0]))
        cell.apply_periodic_bc(step_name, history)
        print("\nRun finite-strain axial load path.")
        print("  Run Abaqus analysis.")
        sys.stdout.flush()
        cell.run_analysis()

        # Keep the ODB open while evaluating every output frame.
        print("  Post-process effective response.")
        sys.stdout.flush()
        cell.open_odb_file()
        try:
            frames = cell.odb.steps[step_name].frames
            curve = []
            for index, frame in enumerate(frames):
                F11 = 1.0 + final_displacement_gradient*frame.frameValue
                stress = cell.average_stress(step_name, frame=index)
                curve.append([F11, stress[0, 0]])
        finally:
            cell.close_odb_file()
        curve = np.asarray(curve)

        print("\n------------------------------")
        print("Finite-strain response ({} frames):".format(len(curve)))
        print("  final F11 = {:.4f}".format(curve[-1, 0]))
        print("  final P11 = {:.6e}".format(curve[-1, 1]))

        result_dir = os.path.join(here, "results")
        if not os.path.isdir(result_dir):
            os.makedirs(result_dir)
        export_to_csv_file(
            curve, os.path.join(result_dir, "response.csv"),
            header="F_11,P_11")

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
        ax.set_title("Finite-strain response of a sheet-Gyroid shell")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(os.path.join(result_dir, "response.png"), dpi=150)
        plt.close(fig)

        print("*** Finished ***")
    finally:
        if gmsh.is_initialized():
            gmsh.finalize()


if __name__ == "__main__":
    # Abaqus replaces stdout; restore the shell that launched the example.
    if sys.__stdout__ is not None:
        sys.stdout = sys.__stdout__
    main()
