# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-or-later

"""Compute the cyclic response of the elastoplastic sandwich cell.

The mixed-dimensional cell combines a soft continuum matrix, two shell
skins, and an embedded BCC beam core. Two fully reversed axial strain cycles
exercise the perfectly plastic structural constituents and produce an effective
hysteresis response.
"""

import os
import sys

import gmsh
import numpy as np
from abaqusConstants import OFF
from sandwich import create_geometry

from abqhom.abqmodel import AbaqusRVE, LoadHistory
from abqhom.gmshmodel import GmshRVE
from abqhom.utils import export_to_csv_file, hill_mandel_residual

here = os.getcwd()


def assign_sections(cell, parameters):
    """Assign the elastic matrix and perfectly plastic constituents."""
    matrix = cell.model.Material(name="MATRIX")
    matrix.Elastic(table=((parameters["matrix_E"],
                           parameters["matrix_nu"]),))

    shell = cell.model.Material(name="SHELL")
    shell.Elastic(table=((parameters["shell_E"],
                          parameters["shell_nu"]),))
    shell.Plastic(table=((parameters["shell_yield"], 0.0),))

    beam = cell.model.Material(name="BEAM")
    beam.Elastic(table=((parameters["beam_E"],
                         parameters["beam_nu"]),))
    beam.Plastic(table=((parameters["beam_yield"], 0.0),))

    cell.assign_solid_section("CONTINUUM_MATRIX", "MATRIX")
    cell.assign_shell_section(
        "SHELL_FACES", "SHELL", parameters["shell_thickness"],
        num_integration_points=5, pre_integrate=False)
    cell.assign_beam_section(
        "BEAM_STRUTS", "BEAM", parameters["beam_radius"],
        pre_integrate=False)


def evaluate_response(cell, step_name, load_path, power_tolerance):
    """Evaluate stresses and Hill--Mandel residuals from an open ODB."""
    response = []
    step = cell.odb.steps[step_name]
    previous_gradient = None
    for frame_index, frame in enumerate(step.frames):
        gradient = load_path(frame.frameValue)
        stress = cell.average_stress(step_name, frame=frame_index)
        scale = max(np.linalg.norm(stress), 1.0e-30)
        asymmetry = np.linalg.norm(stress - stress.T) / scale

        hill_mandel = np.nan
        if previous_gradient is not None:
            # Insert the increment between successive converged frames as the
            # common macroscopic and microscopic test direction.
            increment = gradient - previous_gradient
            macroscopic_work = float(np.sum(stress * increment))
            microscopic_work = cell.average_power_increment(
                step_name, frame=frame_index)
            hill_mandel = hill_mandel_residual(
                macroscopic_work, microscopic_work, power_tolerance)

        response.append(np.concatenate((
            [frame.frameValue, gradient[0, 0]], stress.reshape(-1),
            [asymmetry, hill_mandel])))
        previous_gradient = gradient

    return np.asarray(response)


def evaluate_peeq_history(cell, step_name):
    """Return maximum shell and beam PEEQ values from an open ODB."""
    records = []
    step = cell.odb.steps[step_name]
    instance = cell.odb.rootAssembly.instances[cell.model_name]
    regions = (
        instance.elementSets["ELSET_SHELL_FACES"],
        instance.elementSets["ELSET_BEAM_STRUTS"],
    )

    for frame in step.frames:
        peeq = frame.fieldOutputs["PEEQ"]
        maxima = []
        for region in regions:
            values = peeq.getSubset(region=region).values
            maxima.append(max(float(value.data) for value in values))
        records.append([frame.frameValue] + maxima)

    return np.asarray(records)


def export_results(result_directory, response, peeq_history):
    """Export the response histories and their verification measures."""
    if not os.path.isdir(result_directory):
        os.makedirs(result_directory)

    response_header = (
        "pseudo_time,epsilon_11,sigma_11,sigma_12,sigma_13,"
        "sigma_21,sigma_22,sigma_23,sigma_31,sigma_32,sigma_33,"
        "stress_asymmetry,hill_mandel_residual")
    export_to_csv_file(
        response, os.path.join(result_directory, "cyclic_response.csv"),
        header=response_header)
    export_to_csv_file(
        peeq_history, os.path.join(result_directory, "peeq_history.csv"),
        header="pseudo_time,max_peeq_shell,max_peeq_beam")

    stress_norm = np.linalg.norm(response[:, 2:11], axis=1)
    maximum_stress = np.max(stress_norm)
    if maximum_stress <= 0.0:
        raise RuntimeError("The effective stress is zero in every frame.")
    loaded_frames = stress_norm > 1.0e-6 * maximum_stress
    maximum_asymmetry = np.max(response[loaded_frames, -2])
    maximum_hill_mandel = np.nanmax(response[:, -1])
    maximum_shell_peeq = np.max(peeq_history[:, 1])
    maximum_beam_peeq = np.max(peeq_history[:, 2])

    # Columns: stress asymmetry, Hill--Mandel residual, maximum shell PEEQ,
    # and maximum beam PEEQ.
    metrics = np.asarray([[
        maximum_asymmetry, maximum_hill_mandel,
        maximum_shell_peeq, maximum_beam_peeq,
    ]])
    export_to_csv_file(
        metrics, os.path.join(result_directory, "verification_metrics.csv"),
        header=("maximum_stress_asymmetry,maximum_hill_mandel_residual,"
                "maximum_shell_peeq,maximum_beam_peeq"))

    return {
        "asymmetry": maximum_asymmetry,
        "hill_mandel": maximum_hill_mandel,
        "shell_peeq": maximum_shell_peeq,
        "beam_peeq": maximum_beam_peeq,
    }


def run_case(gmsh_cell, reversal_times, reversal_strains, parameters,
             output_directories):
    """Build, solve, and post-process the cyclic load path."""
    analysis_directory, result_directory, vtk_directory = output_directories
    if not os.path.isdir(analysis_directory):
        os.makedirs(analysis_directory)
    os.chdir(analysis_directory)

    cell = AbaqusRVE(gmsh_cell, custom_types={1: "B31"})
    assign_sections(cell, parameters)

    def load_path(pseudo_time):
        """Evaluate the macroscopic displacement gradient."""
        gradient = np.zeros((3, 3))
        gradient[0, 0] = np.interp(
            pseudo_time, reversal_times, reversal_strains)
        return gradient

    load_history = LoadHistory(reversal_times, load_path)
    time_period = float(reversal_times[-1])
    time_increment = parameters["time_increment"]
    number_of_increments = int(round(time_period / time_increment))
    if not np.isclose(number_of_increments * time_increment, time_period):
        raise ValueError(
            "The fixed increment size must divide the step time exactly.")

    step_name = "LOAD"
    cell.model.StaticStep(
        name=step_name, previous="Initial", nlgeom=OFF,
        timePeriod=time_period, maxNumInc=number_of_increments,
        initialInc=time_increment, minInc=time_increment,
        maxInc=time_increment)

    # Establish periodicity before registering the embedded region. The beam
    # endpoints already share the host nodes at the eight cell corners.
    cell.apply_periodic_bc(step_name, load_history)
    cell.model.EmbeddedRegion(
        "EMBEDDED_BEAM_CORE",
        cell.instance.sets["ELSET_BEAM_STRUTS"],
        cell.instance.sets["ELSET_CONTINUUM_MATRIX"])
    print("\nRun cyclic axial load path.")
    print("  Run Abaqus analysis.")
    sys.stdout.flush()
    cell.run_analysis(field_outputs=["PEEQ"], numCpus=1)

    print("  Post-process effective response and plastic strains.")
    sys.stdout.flush()
    # Reuse one ODB handle for the stress, work, and PEEQ histories.
    cell.open_odb_file()
    try:
        response = evaluate_response(
            cell, step_name, load_path, parameters["power_tolerance"])
        peeq_history = evaluate_peeq_history(cell, step_name)
    finally:
        cell.close_odb_file()
    metrics = export_results(
        result_directory, response, peeq_history)

    print("  Export VTK results.")
    sys.stdout.flush()
    cell.export_to_vtk(
        vtk_directory, field_outputs=["U", "PEEQ"],
        element_sets=["CONTINUUM_MATRIX", "SHELL_FACES", "BEAM_STRUTS"],
        clear=True)

    frame_times = response[:, 0]
    print("\nCyclic effective response:")
    print("  reversal | pseudo-time | epsilon_11 | sigma_11")
    for reversal, (strain, time) in enumerate(
            zip(reversal_strains[1:-1], reversal_times[1:-1]), start=1):
        frame = int(np.argmin(np.abs(frame_times - time)))
        print("  {:8d} | {:11.2f} | {:+.4f}    | {:+.6e}".format(
            reversal, time, strain, response[frame, 2]))
    print("  maximum stress asymmetry     = {:.3e}".format(
        metrics["asymmetry"]))
    print("  maximum Hill--Mandel residual = {:.3e}".format(
        metrics["hill_mandel"]))
    print("  maximum shell PEEQ            = {:.6e}".format(
        metrics["shell_peeq"]))
    print("  maximum beam PEEQ             = {:.6e}".format(
        metrics["beam_peeq"]))


def main():
    """Run the cyclic small-strain plasticity example."""
    original_directory = os.getcwd()
    try:
        # Geometry and material parameters are nondimensional.
        parameters = {
            "matrix_E": 1.0,
            "matrix_nu": 0.3,
            "shell_E": 100.0,
            "shell_nu": 0.3,
            "shell_yield": 0.4,
            "beam_E": 100.0,
            "beam_nu": 0.3,
            "beam_yield": 0.6,
            "shell_thickness": 0.02,
            "beam_radius": 0.05,
            "time_increment": 0.002,
            "power_tolerance": 1.0e-8,
        }
        mesh_parameters = {
            "edge_length": 1.0,
            "matrix_elements_per_edge": 8,
            "elements_per_strut": 16,
            "element_order": 1,
            "algorithm_2d": 6,
            "algorithm_3d": 1,
        }

        # Mesh and loading settings. The reversal times give a constant
        # strain-rate magnitude over the complete piecewise-linear history.
        reversal_strains = np.asarray([
            0.0, 0.03, -0.03, 0.03, -0.03, 0.03, 0.0])
        reversal_times = np.asarray([
            0.0, 0.10, 0.30, 0.50, 0.70, 0.90, 1.00])

        output_directories = (
            os.path.join(here, "abq"),
            os.path.join(here, "results"),
            os.path.join(here, "vtk_output"),
        )

        # ---------------------------------------------------------------
        # Build the mixed-dimensional periodic mesh once.
        model_name = "RVE"
        create_geometry(
            model_name, mesh_parameters["edge_length"],
            mesh_parameters["elements_per_strut"])
        gmsh_cell = GmshRVE(model_name)
        gmsh_cell.add_periodic_mesh_constraints()
        gmsh_cell.generate_mesh(
            element_order=mesh_parameters["element_order"],
            algorithm_2d=mesh_parameters["algorithm_2d"],
            algorithm_3d=mesh_parameters["algorithm_3d"],
            element_size=(mesh_parameters["edge_length"]
                          / mesh_parameters["matrix_elements_per_edge"]))
        gmsh_cell.build_node_sets()
        gmsh_cell.build_element_sets()

        # ---------------------------------------------------------------
        # Solve and post-process the cyclic load path.
        run_case(
            gmsh_cell, reversal_times, reversal_strains, parameters,
            output_directories)

        print("*** Finished cyclic sandwich example ***")
    finally:
        os.chdir(original_directory)
        if gmsh.is_initialized():
            gmsh.finalize()


if __name__ == "__main__":
    # Abaqus replaces stdout. Restore the shell that launched the example.
    if sys.__stdout__ is not None:
        sys.stdout = sys.__stdout__
    main()
