# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-or-later

"""Homogenize a two-phase spherical-inclusion cell at finite strain."""

import os
import sys

import gmsh
import numpy as np
from abaqusConstants import ISOTROPIC, NEO_HOOKE, OFF, ON, VOLUMETRIC_DATA

from abqhom.abqmodel import AbaqusRVE, LoadHistory
from abqhom.gmshmodel import GmshRVE
from abqhom.utils import export_to_csv_file

here = os.getcwd()


def create_geometry(model_name, dx=1.0, dy=1.0, dz=1.0,
                    inclusion_fraction=0.2):
    """Create a cuboidal matrix with a centered spherical inclusion."""
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
    matrix, _ = gmsh.model.occ.cut(
        [(3, box)], [(3, inclusion)], removeTool=False)
    gmsh.model.occ.synchronize()

    group = gmsh.model.add_physical_group(3, [tag for _, tag in matrix])
    gmsh.model.set_physical_name(3, group, "MATRIX")
    group = gmsh.model.add_physical_group(3, [inclusion])
    gmsh.model.set_physical_name(3, group, "INCLUSION")


def main():
    """Run the example workflow."""
    try:
        # ---------------------------------------------------------------
        # Working directory (created next to this script)
        workdir = os.path.join(here, "abq")
        if not os.path.isdir(workdir):
            os.makedirs(workdir)
        os.chdir(workdir)

        # RVE and material parameters (E, nu)
        dx = dy = dz = 10.0
        inclusion_fraction = 0.3       # volume fraction of the inclusion
        parameters = {"MATRIX":    [3.5, 0.35],    # soft matrix
                      "INCLUSION": [230.0, 0.2]}   # stiff inclusion
        final_displacement_gradient = 0.2

        # ---------------------------------------------------------------
        # Build the geometry and impose a matching periodic mesh.
        create_geometry("RVE", dx=dx, dy=dy, dz=dz,
                        inclusion_fraction=inclusion_fraction)
        unit_cell = GmshRVE("RVE")
        unit_cell.add_periodic_mesh_constraints()
        unit_cell.generate_mesh(element_size=2.0)

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
        unit_cell.model.StaticStep(
            name=step_name, previous="Initial", nlgeom=ON, maxNumInc=200,
            initialInc=increment, minInc=increment, maxInc=increment)

        history = LoadHistory(
            np.array([0.0, 1.0]),
            lambda pseudo_time: np.diag(
                [final_displacement_gradient*pseudo_time, 0.0, 0.0]))
        unit_cell.apply_periodic_bc(step_name, history)
        print("\nRun finite-strain axial load path.")
        print("  Run Abaqus analysis.")
        sys.stdout.flush()
        unit_cell.run_analysis(field_outputs=["U", "S"])

        # Evaluate the effective first Piola--Kirchhoff stress frame by frame.
        print("  Post-process effective response.")
        sys.stdout.flush()
        unit_cell.open_odb_file()
        try:
            frames = unit_cell.odb.steps[step_name].frames
            curve = []
            for index, frame in enumerate(frames):
                F11 = 1.0 + final_displacement_gradient*frame.frameValue
                stress = unit_cell.average_stress(step_name, frame=index)
                curve.append([F11, stress[0, 0]])
        finally:
            unit_cell.close_odb_file()
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
        ax.set_title("Finite-strain response of a spherical-inclusion cell")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(os.path.join(result_dir, "response.png"), dpi=150)
        plt.close(fig)

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
