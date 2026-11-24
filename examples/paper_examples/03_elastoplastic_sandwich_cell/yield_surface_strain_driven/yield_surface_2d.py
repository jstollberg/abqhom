# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-or-later

"""Sample sections of the initial and evolved effective yield surfaces.

Every direction uses a fresh Abaqus model and a prescribed macroscopic-strain
history. Initial yield is probed from the stress-free state. Evolved yield is
probed after monotonic preloading and elastic unloading to a prescribed base
strain. Both cases use the same additional plastic-work-density threshold and
PEEQ guard relative to their respective base states.

Run ``abaqus cae noGUI=yield_surface_2d.py`` after selecting the desired
surfaces and planes in the final call to ``main``. Partially completed plane
sections are resumed from their CSV files.
"""

import csv
import os
import sys

import gmsh
import numpy as np

from abqhom.utils import export_to_csv_file
from yield_surface_common import (CONSTITUENT_PARAMETERS, ELASTIC_STRAIN,
                                  MAXIMUM_STRAIN_INCREMENT,
                                  MESH_PARAMETERS,
                                  PLASTIC_WORK_DENSITY_THRESHOLD,
                                  STRAIN_HEADER, STRESS_HEADER,
                                  VOIGT_COMPONENTS, build_gmsh_cell,
                                  elastic_calibration, open_csv_file,
                                  sample_direction, state_definition,
                                  validate_result_configuration)

script_directory = os.getcwd()


# ---------------------------------------------------------------
# Section result persistence

def write_yield_results(path, records):
    """Write sampled yield points and their common base state."""
    header = (["theta", "yield-time", "strain-limit",
               "additional-plastic-work-density"]
              + ["epsilon-{}".format(name)
                 for name in VOIGT_COMPONENTS]
              + ["sigma-{}".format(name)
                 for name in VOIGT_COMPONENTS]
              + ["base-epsilon-{}".format(name)
                 for name in VOIGT_COMPONENTS]
              + ["base-sigma-{}".format(name)
                 for name in VOIGT_COMPONENTS]
              + ["max-unload-delta-peeq-shell",
                 "max-unload-delta-peeq-beam",
                 "max-delta-peeq-shell", "max-delta-peeq-beam",
                 "first-yielding-constituent"])
    with open_csv_file(path, "w") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(header)
        for record in records:
            writer.writerow(
                [record["theta"], record["time"],
                 record["strain_limit"],
                 record["additional_plastic_work_density"]]
                + list(record["strain"])
                + list(record["stress"])
                + list(record["base_strain"])
                + list(record["base_stress"])
                + [record["unload_peeq_shell"],
                   record["unload_peeq_beam"],
                   record["peeq_shell"], record["peeq_beam"],
                   record["first_yielding_constituent"]])


def read_yield_results(path):
    """Read completed yield directions for a resumable sampling run."""
    if not os.path.isfile(path):
        return []

    records = []
    with open_csv_file(path, "r") as csv_file:
        for row in csv.DictReader(csv_file):
            records.append({
                "theta": float(row["theta"]),
                "time": float(row["yield-time"]),
                "strain_limit": float(row["strain-limit"]),
                "additional_plastic_work_density": float(
                    row["additional-plastic-work-density"]),
                "strain": np.asarray([
                    float(row["epsilon-" + name])
                    for name in VOIGT_COMPONENTS]),
                "stress": np.asarray([
                    float(row["sigma-" + name])
                    for name in VOIGT_COMPONENTS]),
                "base_strain": np.asarray([
                    float(row["base-epsilon-" + name])
                    for name in VOIGT_COMPONENTS]),
                "base_stress": np.asarray([
                    float(row["base-sigma-" + name])
                    for name in VOIGT_COMPONENTS]),
                "unload_peeq_shell": float(
                    row["max-unload-delta-peeq-shell"]),
                "unload_peeq_beam": float(
                    row["max-unload-delta-peeq-beam"]),
                "peeq_shell": float(row["max-delta-peeq-shell"]),
                "peeq_beam": float(row["max-delta-peeq-beam"]),
                "first_yielding_constituent": row[
                    "first-yielding-constituent"],
            })
    return records


# ---------------------------------------------------------------
# Plane sampling

def sample_plane(gmsh_cell, parameters, compliance, state_name,
                 preload_strain, base_strain, plane, angles,
                 strain_limits, plastic_work_density_threshold,
                 maximum_strain_increment, result_directory):
    """Sample one resumable normal-stress-plane section."""
    plane_name, first, second = plane
    result_path = os.path.join(result_directory, plane_name + ".csv")

    # A completed prefix of angular directions can be resumed safely only when
    # it belongs to the same ordered angular grid.
    records = read_yield_results(result_path)
    completed_angles = [record["theta"] for record in records]
    if (len(records) > len(angles)
            or not np.allclose(completed_angles, angles[:len(records)])):
        print("Existing {} results do not match the angular grid. "
              "Starting this plane again.".format(plane_name))
        records = []

    print("\nSampling {} {} with {} fixed-increment directions.".format(
        state_name, plane_name, len(angles)))
    if records:
        print("  resuming after {} completed directions".format(
            len(records)))

    for direction_index in range(len(records), len(angles)):
        angle = angles[direction_index]

        # Prescribe a target stress-space direction and use the effective
        # elastic compliance to obtain its proportional strain direction.
        stress_direction = np.zeros(6)
        stress_direction[first] = np.cos(angle)
        stress_direction[second] = np.sin(angle)
        strain_direction = compliance.dot(stress_direction)
        strain_direction /= np.max(np.abs(strain_direction))

        print("  direction {:02d}/{:02d}: theta = {:7.2f} deg".format(
            direction_index + 1, len(angles), np.degrees(angle)))
        sys.stdout.flush()

        # A single clean retry handles occasional solver or file-access
        # failures without silently continuing after a repeated failure.
        for attempt in range(2):
            try:
                state = sample_direction(
                    gmsh_cell, parameters, preload_strain, base_strain,
                    strain_direction, strain_limits,
                    plastic_work_density_threshold,
                    maximum_strain_increment)
                break
            except Exception as error:
                if attempt == 1:
                    raise
                print("    analysis failed once. Rebuilding and retrying: "
                      "{}".format(error))
                sys.stdout.flush()
        state["theta"] = angle
        records.append(state)

        # Persist every completed direction so long sampling runs are
        # resumable.
        write_yield_results(result_path, records)

    print("Completed {}.".format(plane_name))
    sys.stdout.flush()


# ---------------------------------------------------------------
# Script orchestration

def main(selected_state="both", selected_plane=None,
         calibration_only=False):
    """Compute initial and evolved effective yield-surface sections."""
    # Analysis configuration.
    number_of_directions = 36

    original_directory = os.getcwd()
    try:
        gmsh_cell = build_gmsh_cell(MESH_PARAMETERS)

        # The elastic calibration is shared by the initial and evolved sweeps.
        calibration_directory = os.path.join(
            script_directory, "2d_section_results_initial")
        calibration_analysis = os.path.join(
            calibration_directory, "abq_calibration")
        for directory in (calibration_analysis, calibration_directory):
            if not os.path.isdir(directory):
                os.makedirs(directory)
        calibration_results = (
            os.path.join(calibration_directory, "effective_stiffness.csv"),
            os.path.join(calibration_directory, "effective_compliance.csv"),
        )
        validate_result_configuration(
            os.path.join(
                calibration_directory, "calibration_configuration.csv"),
            {
                "analysis": "elastic_calibration",
                "elastic_strain": ELASTIC_STRAIN,
                "materials": CONSTITUENT_PARAMETERS,
                "mesh": MESH_PARAMETERS,
            },
            calibration_results)
        os.chdir(calibration_analysis)
        compliance = elastic_calibration(
            gmsh_cell, CONSTITUENT_PARAMETERS, ELASTIC_STRAIN,
            calibration_directory)
        if calibration_only:
            print("*** Elastic calibration finished ***")
            return

        # The two requested normal-stress planes use the same complete angular
        # grid without invoking symmetry of the cell or effective response.
        if number_of_directions < 4:
            raise ValueError(
                "At least four yield-surface directions are required.")
        angles = np.linspace(
            0.0, 2.0 * np.pi, number_of_directions, endpoint=False)
        stress_planes = (
            ("sigma11_sigma22", 0, 1),
            ("sigma11_sigma33", 0, 2),
        )
        plane_names = tuple(plane[0] for plane in stress_planes)
        if selected_plane is not None and selected_plane not in plane_names:
            raise ValueError(
                "The selected plane must be one of {} or None."
                .format(plane_names))

        selected_state = selected_state.lower()
        if selected_state not in ("initial", "evolved", "both"):
            raise ValueError(
                "The yield state must be initial, evolved, or both.")

        # Each surface has its own Abaqus working directory and resumable CSVs.
        for state_name in ("initial", "evolved"):
            if selected_state not in (state_name, "both"):
                continue
            definition = state_definition(state_name)
            result_directory = os.path.join(
                script_directory,
                "2d_section_results_{}".format(state_name))
            analysis_name = "abq"
            if selected_plane is not None:
                analysis_name += "_" + selected_plane
            analysis_directory = os.path.join(
                result_directory, analysis_name)
            for directory in (analysis_directory, result_directory):
                if not os.path.isdir(directory):
                    os.makedirs(directory)
            result_paths = tuple(
                os.path.join(result_directory, name + ".csv")
                for name in plane_names)
            result_paths += (
                os.path.join(result_directory, "base_strain.csv"),
                os.path.join(result_directory, "base_stress.csv"),
            )
            validate_result_configuration(
                os.path.join(
                    result_directory, "sampling_configuration.csv"),
                {
                    "analysis": "yield_surface_2d",
                    "directions": number_of_directions,
                    "elastic_strain": ELASTIC_STRAIN,
                    "materials": CONSTITUENT_PARAMETERS,
                    "maximum_strain_increment": MAXIMUM_STRAIN_INCREMENT,
                    "mesh": MESH_PARAMETERS,
                    "planes": plane_names,
                    "preload_strain": definition["preload_strain"],
                    "plastic_work_density_threshold": (
                        PLASTIC_WORK_DENSITY_THRESHOLD),
                    "base_strain": definition["base_strain"],
                    "state": state_name,
                    "strain_limits": definition["strain_limits"],
                },
                result_paths)
            os.chdir(analysis_directory)

            print("\n{} yield surface".format(state_name.capitalize()))
            print("  preload strain: {}".format(
                definition["preload_strain"]))
            print("  base strain: {}".format(definition["base_strain"]))
            for plane in stress_planes:
                if (selected_plane is not None
                        and plane[0] != selected_plane):
                    continue
                sample_plane(
                    gmsh_cell, CONSTITUENT_PARAMETERS, compliance,
                    state_name, definition["preload_strain"],
                    definition["base_strain"], plane, angles,
                    definition["strain_limits"],
                    PLASTIC_WORK_DENSITY_THRESHOLD,
                    MAXIMUM_STRAIN_INCREMENT,
                    result_directory=result_directory)

            # Independently reproduced base states must be identical across all
            # probing directions and both sampled planes.
            records = []
            for plane in stress_planes:
                path = os.path.join(
                    result_directory, plane[0] + ".csv")
                if os.path.isfile(path):
                    records.extend(read_yield_results(path))
            if records:
                base_strains = np.asarray([
                    record["base_strain"] for record in records])
                base_stresses = np.asarray([
                    record["base_stress"] for record in records])
                if not np.allclose(
                        base_strains, base_strains[0],
                        rtol=1.0e-6, atol=1.0e-9):
                    raise RuntimeError(
                        "The independently reproduced base states do not "
                        "have consistent effective strains.")
                if not np.allclose(
                        base_stresses, base_stresses[0],
                        rtol=1.0e-6, atol=1.0e-9):
                    raise RuntimeError(
                        "The independently reproduced base states do not "
                        "have consistent effective stresses.")
                export_to_csv_file(
                    base_strains[:1],
                    os.path.join(result_directory, "base_strain.csv"),
                    header=STRAIN_HEADER)
                export_to_csv_file(
                    base_stresses[:1],
                    os.path.join(result_directory, "base_stress.csv"),
                    header=STRESS_HEADER)

        print("*** Yield-surface sampling finished ***")
    finally:
        os.chdir(original_directory)
        if gmsh.is_initialized():
            gmsh.finalize()


if __name__ == "__main__":
    # Abaqus replaces stdout. Restore the shell that launched the example.
    if sys.__stdout__ is not None:
        sys.stdout = sys.__stdout__
    main(selected_state="both", selected_plane=None,
         calibration_only=False)
