# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-only

"""Homogenize a perfectly plastic beam lattice under cyclic loading.

The small-strain example isolates path-dependent material behavior from
geometric nonlinearity and microscopic hardening. Two fully reversed cycles
are applied to a simple-cubic lattice, and the effective hysteresis response
is written to CSV and plotted.
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
from abqhom.utils import export_to_csv_file, voigt_to_tensor

from abaqusConstants import OFF


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
    coords = [(x, y, z) for x in (-h, h) for y in (-h, h) for z in (-h, h)]
    corners = [g.add_point(*c) for c in coords]

    struts = []
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
    gmsh.model.set_physical_name(1, group, "CUBIC")


def main():
    """Run the example workflow."""
    try:
        # ---------------------------------------------------------------
        # Working directory (created next to this script)
        workdir = os.path.join(here, "abq")
        if not os.path.isdir(workdir):
            os.makedirs(workdir)
        os.chdir(workdir)

        # Geometry and perfectly plastic material
        length = 10.0
        radius = 0.5
        element_size = 1.0
        E, nu = 230.0, 0.3
        yield_stress = 1.0

        # Two constant-amplitude cycles at constant absolute strain rate.
        strain_amplitude = 0.01
        reversals = strain_amplitude*np.array(
            [0.0, 1.0, -1.0, 1.0, -1.0, 1.0, 0.0])
        strain_jumps = np.abs(np.diff(reversals))
        reversal_times = np.r_[
            0.0, np.cumsum(strain_jumps)/np.sum(strain_jumps)]


        def load_path(time):
            """Return the prescribed macroscopic strain at ``time``."""
            strain = np.zeros(6)
            strain[0] = np.interp(time, reversal_times, reversals)
            return voigt_to_tensor(strain, "STRAIN")


        history = LoadHistory(reversal_times, load_path)


        # ---------------------------------------------------------------
        # Create geometry and periodic mesh in Gmsh
        model_name = "RUC"
        create_geometry(model_name, length)
        unit_cell = GmshRVE(model_name)
        unit_cell.add_periodic_mesh_constraints()
        unit_cell.generate_mesh(element_order=1, element_size=element_size)

        # Load the model into Abaqus (2-node linear Timoshenko beams)
        unit_cell = AbaqusRVE(unit_cell, custom_types={1: "B31"})

        # Perfectly plastic material
        material = unit_cell.model.Material(name="CUBIC")
        material.Elastic(table=[[E, nu]])
        material.Plastic(table=[[yield_stress, 0.0]])
        unit_cell.assign_beam_section("CUBIC", "CUBIC", radius)

        # Use fixed pseudo-time increments to resolve every linear path segment.
        step_name = "LOAD"
        increment = 0.005
        unit_cell.model.StaticStep(
            name=step_name, previous="Initial", nlgeom=OFF, maxNumInc=1000,
            initialInc=increment, minInc=increment, maxInc=increment)

        # Apply periodic boundary conditions and run the analysis
        unit_cell.apply_periodic_BC(step_name, history)
        unit_cell.run_analysis(field_outputs=["U", "S", "PEEQ"])

        # ---------------------------------------------------------------
        # Effective axial stress--strain response evaluated frame by frame
        unit_cell.open_odb_file()
        try:
            frame_times = [frame.frameValue
                           for frame in unit_cell.odb.steps[step_name].frames]
        finally:
            unit_cell.close_odb_file()

        curve = []
        for index, frame_time in enumerate(frame_times):
            strain = load_path(frame_time)[0, 0]
            stress = unit_cell.average_stress(step_name, frame=index)
            curve.append([strain, stress[0, 0]])
        curve = np.array(curve)

        print("\n------------------------------")
        print("Cyclic effective axial response ({} frames):".format(len(curve)))
        print("  strain amplitude       = {:.4f}".format(strain_amplitude))
        print("  maximum axial stress   = {:.6e}".format(np.max(curve[:, 1])))
        print("  minimum axial stress   = {:.6e}".format(np.min(curve[:, 1])))
        print("  final axial stress     = {:.6e}".format(curve[-1, 1]))

        # Store the response curve and export the deformed lattice to VTK
        result_dir = os.path.join(here, "results")
        if not os.path.isdir(result_dir):
            os.makedirs(result_dir)
        export_to_csv_file(curve, os.path.join(result_dir, "cyclic_response.csv"))

        # ---------------------------------------------------------------
        # Plot the cyclic stress/strain response (matplotlib ships with Abaqus).
        # Use the non-interactive Agg backend so it works under 'noGUI'.
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(6.0, 5.0))
        ax.plot(curve[:, 0], curve[:, 1], "-", color="#1f6feb", lw=1.4)
        ax.axhline(0.0, color="0.6", lw=0.6)
        ax.axvline(0.0, color="0.6", lw=0.6)
        ax.set_xlabel(r"Macroscopic strain  $\bar{\varepsilon}_{11}$")
        ax.set_ylabel(r"Effective stress  $\bar{\sigma}_{11}$")
        ax.set_title("Cyclic response of a perfectly plastic lattice")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(os.path.join(result_dir, "cyclic_response.png"), dpi=150)
        plt.close(fig)

        if PARAQUS_INSTALLED:
            unit_cell.export_to_vtk(os.path.join(here, "vtk_output"),
                                    field_outputs=["U", ("S", "MISES"), "PEEQ"],
                                    clear=True)
        else:
            print("Paraqus is not installed; VTK export was skipped.")

        print("*** Finished ***")
    finally:
        if gmsh.is_initialized():
            gmsh.finalize()


if __name__ == "__main__":
    main()
