# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-or-later

"""Shared helpers for strain-driven effective yield-surface sampling."""

import csv
import errno
import json
import os
import sys
import tempfile
import time

# The entry-point scripts run from this subdirectory. Add the parent example
# directory so they can reuse the common sandwich-cell geometry definition.
example_directory = os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))
if example_directory not in sys.path:
    sys.path.insert(0, example_directory)

import numpy as np
from abaqusConstants import OFF

from abqhom.abqmodel import AbaqusRVE, LoadHistory
from abqhom.gmshmodel import GmshRVE
from abqhom.utils import (export_to_csv_file, read_csv_file,
                          tensor_to_voigt, voigt_to_tensor)
from sandwich import create_geometry

VOIGT_COMPONENTS = ("11", "22", "33", "23", "13", "12")
STRAIN_HEADER = ",".join(
    "epsilon_" + component for component in VOIGT_COMPONENTS)
STRESS_HEADER = ",".join(
    "sigma_" + component for component in VOIGT_COMPONENTS)

# These settings define the same model and operational yield criterion for
# the two-dimensional sections and the three-dimensional surface.
DISSIPATION_RATIO_THRESHOLD = 1.0e-4
MAXIMUM_STRAIN_INCREMENT = 6.0e-5
ELASTIC_STRAIN = 1.0e-5
CONSTITUENT_PARAMETERS = {
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
}
MESH_PARAMETERS = {
    "edge_length": 1.0,
    "matrix_elements_per_edge": 8,
    "elements_per_strut": 16,
    "element_order": 1,
    "algorithm_2d": 6,
    "algorithm_3d": 1,
}


def open_csv_file(path, mode):
    """Open a CSV file correctly under Python 2.7 and Python 3."""
    if sys.version_info[0] < 3:
        return open(path, mode + "b")
    return open(path, mode, newline="")


def _open_csv_descriptor(descriptor):
    """Open a temporary CSV descriptor for the active Python version."""
    if sys.version_info[0] < 3:
        return os.fdopen(descriptor, "wb")
    return os.fdopen(descriptor, "w", newline="")


# ---------------------------------------------------------------
# Constituent definitions and effective elastic calibration

def _assign_sections(cell, parameters):
    """Assign the elastic matrix and perfectly plastic constituents."""
    # Plasticity is confined to the shell and beam constituents. The soft
    # matrix remains elastic throughout every preparation and probing path.
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

    # Active section integration resolves yielding through the shell
    # thickness and across the circular beam section.
    cell.assign_shell_section(
        "SHELL_FACES", "SHELL", parameters["shell_thickness"],
        num_integration_points=5, pre_integrate=False)
    cell.assign_beam_section(
        "BEAM_STRUTS", "BEAM", parameters["beam_radius"],
        pre_integrate=False)


def _compute_elastic_stiffness(gmsh_cell, parameters, strain_magnitude):
    """Compute the complete effective elastic stiffness in Voigt notation."""
    stiffness = np.zeros((6, 6))
    component_names = ("xx", "yy", "zz", "yz", "xz", "xy")

    # Each independent macroscopic strain case supplies one stiffness column.
    for component, component_name in enumerate(component_names):
        print("Run elastic strain case {} ({}).".format(
            component + 1, component_name))
        sys.stdout.flush()

        def strain_path(pseudo_time, index=component):
            """Return one ramped macroscopic strain component."""
            strain = np.zeros(6)
            strain[index] = pseudo_time * strain_magnitude
            return voigt_to_tensor(strain, "STRAIN")

        history = LoadHistory(np.asarray([0.0, 1.0]), strain_path)
        cell = AbaqusRVE(gmsh_cell, custom_types={1: "B31"})
        _assign_sections(cell, parameters)
        step_name = "ELASTIC_LOAD"
        cell.model.StaticStep(
            name=step_name, previous="Initial", nlgeom=OFF,
            timePeriod=1.0, maxNumInc=1,
            initialInc=1.0, minInc=1.0, maxInc=1.0)
        cell.apply_periodic_bc(step_name, history)
        cell.model.EmbeddedRegion(
            "EMBEDDED_BEAM_CORE",
            cell.instance.sets["ELSET_BEAM_STRUTS"],
            cell.instance.sets["ELSET_CONTINUUM_MATRIX"])
        print("  Run Abaqus analysis.")
        sys.stdout.flush()
        cell.run_analysis(numCpus=1)
        stiffness[:, component] = (
            tensor_to_voigt(cell.average_stress(step_name), "STRESS")
            / strain_magnitude)

    return stiffness


def elastic_calibration(gmsh_cell, parameters, strain_magnitude,
                        result_directory):
    """Return the effective elastic compliance used for stress directions."""
    stiffness_path = os.path.join(
        result_directory, "effective_stiffness.csv")
    compliance_path = os.path.join(
        result_directory, "effective_compliance.csv")

    # Reuse the six elastic calibration analyses across both yield surfaces.
    if os.path.isfile(stiffness_path) and os.path.isfile(compliance_path):
        stiffness = read_csv_file(stiffness_path)
        compliance = read_csv_file(compliance_path)
        print("Reusing the existing effective elastic calibration.")
    else:
        stiffness = _compute_elastic_stiffness(
            gmsh_cell, parameters, strain_magnitude)
        compliance = None

    if stiffness.shape != (6, 6) or not np.all(np.isfinite(stiffness)):
        raise RuntimeError(
            "The effective stiffness must be a finite 6 x 6 matrix.")

    # The symmetric part defines the stress-to-strain direction map. The raw
    # matrix is retained in the output so its numerical symmetry can be
    # checked.
    stiffness_symmetric = 0.5 * (stiffness + stiffness.T)
    symmetry_residual = (np.linalg.norm(stiffness - stiffness.T)
                         / np.linalg.norm(stiffness))
    eigenvalues = np.linalg.eigvalsh(stiffness_symmetric)
    if np.min(eigenvalues) <= 0.0:
        raise RuntimeError("The effective elastic stiffness is not positive.")
    if compliance is None:
        compliance = np.linalg.inv(stiffness_symmetric)
        export_to_csv_file(stiffness, stiffness_path)
        export_to_csv_file(compliance, compliance_path)
    if compliance.shape != (6, 6) or not np.all(np.isfinite(compliance)):
        raise RuntimeError(
            "The effective compliance must be a finite 6 x 6 matrix.")
    if not np.allclose(
            stiffness_symmetric.dot(compliance), np.eye(6),
            rtol=1.0e-6, atol=1.0e-9):
        raise RuntimeError(
            "The cached effective stiffness and compliance are inconsistent.")
    print("Elastic stiffness symmetry residual: {:.3e}".format(
        symmetry_residual))
    print("Elastic stiffness eigenvalues: {}".format(eigenvalues))
    sys.stdout.flush()
    return compliance


# ---------------------------------------------------------------
# Incremental yield detection from the converged Abaqus output

def _matching_frame_index(frame_times, pseudo_time):
    """Return the field frame matching a history-output time."""
    index = int(np.argmin(np.abs(frame_times - pseudo_time)))
    if not np.isclose(
            frame_times[index], pseudo_time, rtol=1.0e-10,
            atol=1.0e-12):
        raise RuntimeError(
            "Field and history output do not share the required increment "
            "at pseudo-time {}.".format(pseudo_time))
    return index


def _constituent_peeq_increment_data(cell, step_name, lower_time, upper_time,
                                     tolerance=1.0e-12):
    """Return constituent PEEQ increments and first-yield times.

    PEEQ increments are measured from the first frame of the requested step.
    Because PEEQ is non-decreasing, the first frame whose increment exceeds
    the tolerance is located by a binary search rather than by reading every
    output frame.
    """
    close_odb = cell.odb is None
    cell.open_odb_file()
    try:
        step = cell.odb.steps[step_name]
        instance = cell.odb.rootAssembly.instances[cell.model_name]
        constituent_sets = (
            instance.elementSets["ELSET_SHELL_FACES"],
            instance.elementSets["ELSET_BEAM_STRUTS"],
        )

        # The interpolation below assumes that field and history output share
        # the converged increments bracketing the work-density threshold.
        frames = step.frames
        frame_times = np.asarray(
            [frame.frameValue for frame in frames], dtype=float)
        lower_frame = _matching_frame_index(frame_times, lower_time)
        upper_frame = _matching_frame_index(frame_times, upper_time)

        # Subtracting the first frame makes the same operation valid for the
        # stress-free initial state and the plastically preloaded base state.
        baseline_output = frames[0].fieldOutputs["PEEQ"]
        baseline = tuple(
            np.asarray([
                float(value.data) for value in
                baseline_output.getSubset(region=region).values])
            for region in constituent_sets)

        # PEEQ extraction is cached because the binary searches revisit frames.
        cache = {}

        def constituent_increment(frame_index, constituent):
            """Return the maximum PEEQ increment for one constituent."""
            key = (frame_index, constituent)
            if key not in cache:
                peeq = frames[frame_index].fieldOutputs["PEEQ"]
                values = peeq.getSubset(
                    region=constituent_sets[constituent]).values
                current = np.asarray([
                    float(value.data) for value in values])
                if current.shape != baseline[constituent].shape:
                    raise RuntimeError(
                        "The constituent PEEQ output changed size during "
                        "probing.")
                cache[key] = max(
                    0.0, float(np.max(current - baseline[constituent])))
            return cache[key]

        # Locate the first frame with additional plastic straining separately
        # for shells and beams to identify the first-yielding constituent.
        onset_times = [np.nan, np.nan]
        for constituent in range(len(constituent_sets)):
            if (constituent_increment(len(frames) - 1, constituent)
                    <= tolerance):
                continue
            left = 0
            right = len(frames) - 1
            while left < right:
                middle = (left + right) // 2
                if constituent_increment(middle, constituent) > tolerance:
                    right = middle
                else:
                    left = middle + 1
            onset_times[constituent] = frame_times[left]

        peeq = np.asarray([
            [constituent_increment(lower_frame, constituent)
             for constituent in range(len(constituent_sets))],
            [constituent_increment(upper_frame, constituent)
             for constituent in range(len(constituent_sets))],
        ])
    finally:
        if close_odb:
            cell.close_odb_file()
    return peeq, np.asarray(onset_times)


def _yield_state(cell, step_name, base_strain, strain_direction,
                 strain_limit, dissipation_ratio_threshold):
    """Interpolate yield from work increments accumulated during one probe."""
    close_odb = cell.odb is None
    cell.open_odb_file()
    try:
        step = cell.odb.steps[step_name]
        history_region = step.historyRegions["Assembly ASSEMBLY"]
        plastic = np.asarray(
            history_region.historyOutputs["ALLPD"].data, dtype=float)
        internal = np.asarray(
            history_region.historyOutputs["ALLIE"].data, dtype=float)
        if not np.allclose(plastic[:, 0], internal[:, 0]):
            raise RuntimeError(
                "ALLPD and ALLIE are not available at identical times.")

        # Only work increments accumulated from the probe base state enter the
        # operational yield criterion. Preload values are subtracted for the
        # evolved surface.
        plastic_increment = np.maximum(
            plastic[:, 1] - plastic[0, 1], 0.0)
        internal_increment = internal[:, 1] - internal[0, 1]
        dissipation_ratio = (
            plastic_increment
            / np.maximum(np.abs(internal_increment), 1.0e-30))

        # First bracket the prescribed dissipation ratio in pseudo-time.
        yielded = np.flatnonzero(
            dissipation_ratio >= dissipation_ratio_threshold)
        if yielded.size == 0:
            return None
        upper = int(yielded[0])
        if upper == 0:
            raise RuntimeError(
                "The initial frame satisfies the yield criterion.")
        lower = upper - 1

        # The PEEQ guard prevents an energy-ratio crossing from being accepted
        # before either structural constituent has accumulated plastic strain.
        _, onset_times = _constituent_peeq_increment_data(
            cell, step_name, plastic[0, 0], plastic[-1, 0])
        if np.all(np.isnan(onset_times)):
            return None
        plastic_onset = np.nanmin(onset_times)
        if plastic[upper, 0] < plastic_onset - 1.0e-12:
            onset_frame = int(np.searchsorted(
                plastic[:, 0], plastic_onset - 1.0e-12))
            eligible = np.flatnonzero(
                (np.arange(len(dissipation_ratio)) >= onset_frame)
                & (dissipation_ratio >= dissipation_ratio_threshold))
            if eligible.size == 0:
                return None
            upper = int(eligible[0])
            lower = upper - 1

        # Interpolate the crossing within the two converged bracket frames.
        denominator = (dissipation_ratio[upper]
                       - dissipation_ratio[lower])
        if (plastic[upper, 0] <= plastic_onset + 1.0e-12
                and dissipation_ratio[lower] >= dissipation_ratio_threshold):
            fraction = 1.0
        else:
            fraction = ((dissipation_ratio_threshold
                         - dissipation_ratio[lower]) / denominator
                        if abs(denominator) > 1.0e-30 else 1.0)
        yield_time = (plastic[lower, 0]
                      + fraction * (plastic[upper, 0] - plastic[lower, 0]))

        # Apply the same interpolation fraction to stress and PEEQ output.
        peeq, _ = _constituent_peeq_increment_data(
            cell, step_name, plastic[lower, 0], plastic[upper, 0])
        frame_times = np.asarray(
            [frame.frameValue for frame in step.frames], dtype=float)
        lower_frame = _matching_frame_index(
            frame_times, plastic[lower, 0])
        upper_frame = _matching_frame_index(
            frame_times, plastic[upper, 0])
        lower_stress = tensor_to_voigt(
            cell.average_stress(step_name, lower_frame), "STRESS")
        upper_stress = tensor_to_voigt(
            cell.average_stress(step_name, upper_frame), "STRESS")
        stress = lower_stress + fraction * (upper_stress - lower_stress)
        strain = base_strain + yield_time * strain_limit * strain_direction
        peeq_yield = peeq[0] + fraction * (peeq[1] - peeq[0])

        # Report which structural constituent first activated the PEEQ guard.
        if np.isnan(onset_times[1]) or (
                not np.isnan(onset_times[0])
                and onset_times[0] < onset_times[1] - 1.0e-12):
            first_yielding_constituent = "SHELL"
        elif np.isnan(onset_times[0]) or (
                onset_times[1] < onset_times[0] - 1.0e-12):
            first_yielding_constituent = "BEAM"
        else:
            first_yielding_constituent = "BOTH"

        return {
            "time": yield_time,
            "strain": strain,
            "stress": stress,
            "dissipation_ratio": dissipation_ratio_threshold,
            "peeq_shell": peeq_yield[0],
            "peeq_beam": peeq_yield[1],
            "first_yielding_constituent": first_yielding_constituent,
        }
    finally:
        if close_odb:
            cell.close_odb_file()


def sample_direction(gmsh_cell, parameters, preload_strain, base_strain_target,
                     strain_direction, strain_limits,
                     dissipation_ratio_threshold,
                     maximum_strain_increment):
    """Prepare one RVE state and probe a proportional strain direction."""
    preload_strain = np.asarray(preload_strain, dtype=float)
    base_strain_target = np.asarray(base_strain_target, dtype=float)
    prepared_state = np.any(np.abs(preload_strain) > 0.0)

    # Increase the probe range only when the preceding range does not reach
    # yield. Every attempt starts from a fresh RVE with an independent history.
    for strain_limit in strain_limits:
        number_of_increments = int(np.ceil(
            strain_limit / maximum_strain_increment))
        probe_increment = 1.0 / number_of_increments

        def strain_path(pseudo_time, direction=strain_direction,
                        limit=strain_limit):
            """Return preloading, unloading, and probing strains."""
            if not prepared_state:
                strain = pseudo_time * limit * direction
            elif pseudo_time <= 1.0:
                strain = pseudo_time * preload_strain
            elif pseudo_time <= 2.0:
                strain = (2.0 - pseudo_time) * preload_strain
                strain += (pseudo_time - 1.0) * base_strain_target
            else:
                strain = base_strain_target.copy()
                strain += (pseudo_time - 2.0) * limit * direction
            return voigt_to_tensor(strain, "STRAIN")

        cell = AbaqusRVE(gmsh_cell, custom_types={1: "B31"})
        _assign_sections(cell, parameters)

        # The initial surface starts directly with the probe. The evolved
        # surface first reproduces its monotonic preload and elastic unloading.
        if prepared_state:
            preload_step = "MONOTONIC_PRELOAD"
            unload_step = "ELASTIC_UNLOAD"
            probe_step = "PROBE"
            cell.model.StaticStep(
                name=preload_step, previous="Initial", nlgeom=OFF,
                timePeriod=1.0, maxNumInc=100,
                initialInc=0.01, minInc=0.01, maxInc=0.01)
            cell.model.StaticStep(
                name=unload_step, previous=preload_step, nlgeom=OFF,
                timePeriod=1.0, maxNumInc=10,
                initialInc=0.1, minInc=0.1, maxInc=0.1)
            cell.model.StaticStep(
                name=probe_step, previous=unload_step, nlgeom=OFF,
                timePeriod=1.0, maxNumInc=number_of_increments,
                initialInc=probe_increment, minInc=probe_increment,
                maxInc=probe_increment)
            history = LoadHistory(
                np.asarray([0.0, 1.0, 2.0, 3.0]), strain_path)
            first_step = preload_step
        else:
            unload_step = None
            probe_step = "PROBE"
            cell.model.StaticStep(
                name=probe_step, previous="Initial", nlgeom=OFF,
                timePeriod=1.0, maxNumInc=number_of_increments,
                initialInc=probe_increment, minInc=probe_increment,
                maxInc=probe_increment)
            history = LoadHistory(np.asarray([0.0, 1.0]), strain_path)
            first_step = probe_step

        # Applying the complete macroscopic-strain history at its first step
        # propagates the periodic displacement conditions through all steps.
        # The conforming beam endpoints share the host nodes at the cell
        # corners, so periodicity precedes embedding of the interior nodes.
        cell.apply_periodic_bc(first_step, history)
        cell.model.EmbeddedRegion(
            "EMBEDDED_BEAM_CORE",
            cell.instance.sets["ELSET_BEAM_STRUTS"],
            cell.instance.sets["ELSET_CONTINUUM_MATRIX"])
        cell.run_analysis(
            field_outputs=["PEEQ"],
            history_outputs=["ALLPD", "ALLIE"], numCpus=1)

        # Reuse one ODB handle for unloading checks and yield extraction.
        cell.open_odb_file()
        try:
            # Confirm that unloading to the effective base strain did not
            # produce additional plastic strain.
            if prepared_state:
                preload_peeq, _ = _constituent_peeq_increment_data(
                    cell, preload_step, 0.0, 1.0)
                if np.max(preload_peeq[-1]) <= 1.0e-12:
                    raise RuntimeError(
                        "The preload did not produce plastic strain.")

                unload_peeq, _ = _constituent_peeq_increment_data(
                    cell, unload_step, 0.0, 1.0)
                maximum_unload_peeq = unload_peeq[-1]
                if np.max(maximum_unload_peeq) > 1.0e-10:
                    raise RuntimeError(
                        "The common base state was not reached by elastic "
                        "unloading. Maximum constituent PEEQ increments are "
                        "{}.".format(maximum_unload_peeq))
                base_stress = tensor_to_voigt(
                    cell.average_stress(unload_step), "STRESS")
                base_strain = base_strain_target.copy()
            else:
                maximum_unload_peeq = np.zeros(2)
                base_stress = np.zeros(6)
                base_strain = np.zeros(6)

            # Extract the first acceptable yield point from the probe.
            state = _yield_state(
                cell, probe_step, base_strain, strain_direction,
                strain_limit, dissipation_ratio_threshold)
        finally:
            cell.close_odb_file()
        if state is not None:
            state.update({
                "strain_limit": strain_limit,
                "base_strain": base_strain,
                "base_stress": base_stress,
                "unload_peeq_shell": maximum_unload_peeq[0],
                "unload_peeq_beam": maximum_unload_peeq[1],
                # The caller may remove this temporary database after the
                # sampled state has been persisted.
                "odb_path": os.path.abspath(cell.job_name + ".odb"),
            })
            return state
    raise RuntimeError(
        "The yield criterion was not reached up to a maximum strain of {}."
        .format(strain_limits[-1]))


# ---------------------------------------------------------------
# Resumable configuration records

def _configuration_value(value):
    """Return a JSON-compatible configuration value."""
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {
            key: _configuration_value(item)
            for key, item in sorted(value.items())
        }
    if isinstance(value, (list, tuple)):
        return [_configuration_value(item) for item in value]
    return value


def _acquire_configuration_lock(path, timeout=30.0):
    """Acquire an interprocess lock for one configuration file."""
    lock_path = path + ".lock"
    deadline = time.time() + timeout
    while True:
        try:
            os.mkdir(lock_path)
            return lock_path
        except OSError as error:
            if error.errno != errno.EEXIST:
                raise
            if time.time() >= deadline:
                raise RuntimeError(
                    "Timed out waiting for the result-configuration lock "
                    "'{}'. Remove a stale lock after confirming that no "
                    "worker is running.".format(lock_path))
            time.sleep(0.05)


def _validate_result_configuration(path, configuration, result_paths):
    """Create or validate one configuration while holding its lock."""
    expected = {
        key: json.dumps(
            _configuration_value(value), sort_keys=True,
            separators=(",", ":"))
        for key, value in sorted(configuration.items())
    }

    if os.path.isfile(path):
        with open_csv_file(path, "r") as csv_file:
            reader = csv.DictReader(csv_file)
            rows = list(reader)
        actual = {row["setting"]: row["value"] for row in rows}
        if len(actual) != len(rows):
            raise RuntimeError(
                "The result configuration '{}' contains duplicate settings."
                .format(path))
        if actual != expected:
            changed = sorted(
                set(actual).symmetric_difference(expected)
                | {key for key in set(actual).intersection(expected)
                   if actual[key] != expected[key]})
            raise RuntimeError(
                "Existing results use a different configuration for {}. "
                "Remove the corresponding generated results before rerunning."
                .format(", ".join(changed)))
        return

    if any(os.path.isfile(result_path) for result_path in result_paths):
        raise RuntimeError(
            "Existing result files have no configuration record. Remove the "
            "corresponding generated results before rerunning.")

    directory = os.path.dirname(os.path.abspath(path))
    descriptor, temporary_path = tempfile.mkstemp(
        prefix="." + os.path.basename(path) + ".", suffix=".tmp",
        dir=directory, text=True)
    try:
        with _open_csv_descriptor(descriptor) as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(["setting", "value"])
            for key in sorted(expected):
                writer.writerow([key, expected[key]])
        os.rename(temporary_path, path)
    finally:
        if os.path.isfile(temporary_path):
            os.remove(temporary_path)


def validate_result_configuration(path, configuration, result_paths=()):
    """Create or validate the configuration of resumable result files."""
    lock_path = _acquire_configuration_lock(path)
    try:
        _validate_result_configuration(path, configuration, result_paths)
    finally:
        os.rmdir(lock_path)


# ---------------------------------------------------------------
# Geometry and state definitions

def build_gmsh_cell(mesh_parameters):
    """Build the periodic mesh shared by the yield-surface analyses."""
    model_name = "RVE"
    edge_length = mesh_parameters["edge_length"]
    create_geometry(
        model_name, edge_length, mesh_parameters["elements_per_strut"])
    gmsh_cell = GmshRVE(model_name)

    # Enforce matching periodic meshes before generating the final node and
    # element sets consumed by the Abaqus model.
    gmsh_cell.add_periodic_mesh_constraints()
    gmsh_cell.generate_mesh(
        element_order=mesh_parameters["element_order"],
        algorithm_2d=mesh_parameters["algorithm_2d"],
        algorithm_3d=mesh_parameters["algorithm_3d"],
        element_size=(edge_length
                      / mesh_parameters["matrix_elements_per_edge"]))
    gmsh_cell.build_node_sets()
    gmsh_cell.build_element_sets()
    return gmsh_cell


def state_definition(state_name):
    """Return the preparation path and probe limits for one yield surface."""
    zero_strain = np.zeros(6)
    definitions = {
        # No preparation is required for the stress-free initial surface.
        "initial": {
            "preload_strain": zero_strain,
            "base_strain": zero_strain,
            "strain_limits": (0.006, 0.012, 0.024, 0.03, 0.06),
        },
        # Prepare the evolved surface by monotonic strain preloading followed
        # by elastic unloading to the specified base strain.
        "evolved": {
            "preload_strain": np.asarray(
                [0.03, 0.0, 0.0,
                 0.0, 0.0, 0.0]),
            "base_strain": np.asarray(
                [0.027, 0.0, 0.0,
                 0.0, 0.0, 0.0]),
            "strain_limits": (0.003, 0.006, 0.012, 0.024),
        },
    }
    if state_name not in definitions:
        raise ValueError("The yield state must be initial or evolved.")
    return definitions[state_name]
