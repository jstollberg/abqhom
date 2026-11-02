# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-only

"""Homogenize a two-phase spherical-inclusion cell at finite strain.

Automatic boundary sampling is unsuitable for the trimmed OpenCASCADE
geometry. The example therefore classifies the cuboidal boundary entities
geometrically and supplies them explicitly to ``GmshRVE`` before computing a
monotonic finite-strain response.
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

from abaqusConstants import ISOTROPIC, NEO_HOOKE, OFF, ON, VOLUMETRIC_DATA

from abqhom.abqmodel import AbaqusRVE, LoadHistory
from abqhom.gmshmodel import Boundary, GmshRVE
from abqhom.utils import export_to_csv_file


def _cuboidal_boundary(dx, dy, dz):
    """Classify the cuboid entities without relying on their Gmsh tags."""
    x_min, x_max = -0.5*dx, 0.5*dx
    y_min, y_max = -0.5*dy, 0.5*dy
    z_min, z_max = -0.5*dz, 0.5*dz
    group_identifiers = {
        "POINT_1": (x_min, y_min, z_min),
        "POINT_2": (x_max, y_min, z_min),
        "POINT_3": (x_max, y_max, z_min),
        "POINT_4": (x_min, y_max, z_min),
        "POINT_5": (x_min, y_min, z_max),
        "POINT_6": (x_max, y_min, z_max),
        "POINT_7": (x_max, y_max, z_max),
        "POINT_8": (x_min, y_max, z_max),
        "EDGE_1":  (None, y_min, z_min),
        "EDGE_2":  (x_max, None, z_min),
        "EDGE_3":  (None, y_max, z_min),
        "EDGE_4":  (x_min, None, z_min),
        "EDGE_5":  (None, y_min, z_max),
        "EDGE_6":  (x_max, None, z_max),
        "EDGE_7":  (None, y_max, z_max),
        "EDGE_8":  (x_min, None, z_max),
        "EDGE_9":  (x_min, y_min, None),
        "EDGE_10": (x_max, y_min, None),
        "EDGE_11": (x_max, y_max, None),
        "EDGE_12": (x_min, y_max, None),
        "FACE_1":  (None, None, z_min),
        "FACE_2":  (None, None, z_max),
        "FACE_3":  (None, y_min, None),
        "FACE_4":  (None, y_max, None),
        "FACE_5":  (x_min, None, None),
        "FACE_6":  (x_max, None, None),
    }

    tolerance = gmsh.option.get_number("Geometry.Tolerance")
    boundary = {}
    for group_name, identifiers in group_identifiers.items():
        dimension = identifiers.count(None)
        boundary[group_name] = []
        for dimtag in gmsh.model.get_entities(dimension):
            bounding_box = gmsh.model.get_bounding_box(*dimtag)
            lower = bounding_box[:3]
            upper = bounding_box[3:]
            matches = all(
                value is None or np.allclose(
                    (lower[index], upper[index]), value,
                    atol=tolerance, rtol=0.0)
                for index, value in enumerate(identifiers))
            if matches:
                boundary[group_name].append(dimtag)
    return boundary


def create_geometry(model_name, dx=1.0, dy=1.0, dz=1.0,
                    inclusion_fraction=0.2):
    """Create a cuboidal matrix with a centered spherical inclusion."""

    # Compute the radius of the inclusion
    radius = ((0.75*inclusion_fraction*dx*dy*dz)/np.pi)**(1.0/3.0)
    if 2*radius >= min(dx, dy, dz):
        raise ValueError("Inclusion is bigger than the RVE.")

    # Initialize Gmsh and add a new model
    if not gmsh.is_initialized():
        gmsh.initialize()
    gmsh.model.add(model_name)

    gmsh.option.set_number("Geometry.Tolerance", 1e-5)

    # Center point should be coordinate system origin
    x = y = z = 0.0

    # Create the model
    box = gmsh.model.occ.add_box(
        x - dx/2.0, y - dy/2.0, z - dz/2.0, dx, dy, dz)
    inclusion = gmsh.model.occ.add_sphere(x, y, z, radius)
    matrix, _ = gmsh.model.occ.cut(
        [(3, box)], [(3, inclusion)], removeTool=False)
    gmsh.model.occ.synchronize()

    # Automatic boundary generation fails for the trimmed spherical surface.
    # Classify the outer cuboid geometrically instead of using entity tags.
    boundary = _cuboidal_boundary(dx, dy, dz)

    # Add physical groups
    group = gmsh.model.add_physical_group(3, [tag for _, tag in matrix])
    gmsh.model.set_physical_name(3, group, "MATRIX")
    group = gmsh.model.add_physical_group(3, [inclusion])
    gmsh.model.set_physical_name(3, group, "INCLUSION")

    return boundary


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

        # RVE and material parameters (E, nu)
        dx = dy = dz = 10.0
        inclusion_fraction = 0.3       # volume fraction of the inclusion
        parameters = {"MATRIX":    [3.5, 0.35],    # soft matrix
                      "INCLUSION": [230.0, 0.2]}   # stiff inclusion
        final_displacement_gradient = 0.2

        # ---------------------------------------------------------------
        # Build the geometry with a manually supplied boundary and a periodic mesh.
        boundary = create_geometry("RVE", dx=dx, dy=dy, dz=dz,
                                   inclusion_fraction=inclusion_fraction)
        unit_cell = GmshRVE("RVE", Boundary("RVE", boundary))
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
        unit_cell.apply_periodic_BC(step_name, history)
        unit_cell.run_analysis(field_outputs=["U", "S"])

        # Evaluate the effective first Piola--Kirchhoff stress frame by frame.
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
        ax.set_ylabel(r"Effective first Piola-Kirchhoff stress  $\bar{P}_{11}$")
        ax.set_title("Finite-strain response of a spherical-inclusion cell")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(os.path.join(result_dir, "response.png"), dpi=150)
        plt.close(fig)

        unit_cell.export_to_vtk(
            os.path.join(here, "vtk_output"),
            field_outputs=["U", ("S", "MISES")],
            element_sets=["MATRIX", "INCLUSION"], clear=True)

        print("*** Finished ***")
    finally:
        if gmsh.is_initialized():
            gmsh.finalize()


if __name__ == "__main__":
    main()
