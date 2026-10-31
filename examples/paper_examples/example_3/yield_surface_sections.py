# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-only

"""Sample sections of the initial and evolved effective yield surfaces.

Every direction uses a fresh RVE. Initial yield is probed from the stress-free
state. Evolved yield is probed after monotonic preloading and elastic unloading
to a common base state. Both cases use the same incremental dissipation
criterion and PEEQ guard relative to their respective base states.

Run ``abaqus cae noGUI=yield_surface_sections.py`` after selecting the desired
surfaces and planes in the configuration block of ``main``. Partially
completed plane sections are resumed from their CSV files.
"""

import csv
import inspect
import json
import os
import sys
import tempfile
import time

import numpy as np

# ---------------------------------------------------------------
# Locate this example and put the shipped abqhom source first on the
# path. Under 'abaqus cae noGUI=yield_surface_sections.py' __file__ is
# undefined, so the script location is obtained via inspect instead.
here = os.path.dirname(os.path.abspath(
    inspect.getframeinfo(inspect.currentframe()).filename))
abqhom_path = os.path.abspath(os.path.join(
    here, os.pardir, os.pardir, os.pardir, "src"))
if abqhom_path not in sys.path:
    sys.path.insert(0, abqhom_path)

import gmsh
from abaqusConstants import OFF

from abqhom.abqmodel import AbaqusRVE, LoadHistory
from abqhom.gmshmodel import GmshRVE
from abqhom.utils import tensor_to_voigt, voigt_to_tensor
from sandwich import create_geometry


# ---------------------------------------------------------------
# Constituent definitions and effective elastic calibration

def assign_sections(cell, parameters):
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


def compute_elastic_stiffness(gmsh_cell, parameters, strain_magnitude):
    """Compute the complete effective elastic stiffness in Voigt notation."""
    stiffness = np.zeros((6, 6))

    # Each independent macroscopic strain case supplies one stiffness column.
    for component in range(6):
        print("Elastic stiffness case {}/6".format(component + 1))
        sys.stdout.flush()

        def strain_path(time, index=component):
            """Return one ramped macroscopic strain component."""
            strain = np.zeros(6)
            strain[index] = time * strain_magnitude
            return voigt_to_tensor(strain, "STRAIN")

        history = LoadHistory(np.asarray([0.0, 1.0]), strain_path)
        cell = AbaqusRVE(gmsh_cell, custom_types={1: "B31"})
        assign_sections(cell, parameters)
        step_name = "ELASTIC_LOAD"
        cell.model.StaticStep(
            name=step_name, previous="Initial", nlgeom=OFF,
            timePeriod=1.0, maxNumInc=1,
            initialInc=1.0, minInc=1.0, maxInc=1.0)
        cell.apply_periodic_BC(step_name, history)
        cell.model.EmbeddedRegion(
            "EMBEDDED_BEAM_CORE",
            cell.instance.sets["ELSET_BEAM_STRUTS"],
            cell.instance.sets["ELSET_CONTINUUM_MATRIX"])
        cell.run_analysis(numCpus=1)
        stiffness[:, component] = (
            tensor_to_voigt(cell.average_stress(step_name), "STRESS")
            / strain_magnitude)

    return stiffness


# ---------------------------------------------------------------
# Incremental yield detection from the converged Abaqus output

def constituent_peeq_increment_data(cell, step_name, lower_time, upper_time,
                                    tolerance=1.0e-12):
    """Return constituent PEEQ increments and first-yield times.

    PEEQ increments are measured from the first frame of the requested step.
    Because PEEQ is non-decreasing, the first frame with a nonzero increment
    is located by a binary search rather than by reading every output frame.
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

        # Frame indices are obtained from the requested pseudo-times because
        # energy histories and field output can use different sampling grids.
        frames = step.frames
        frame_times = np.asarray(
            [frame.frameValue for frame in frames], dtype=float)
        lower_frame = int(np.argmin(np.abs(frame_times - lower_time)))
        upper_frame = int(np.argmin(np.abs(frame_times - upper_time)))

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
        # for shells and beams to identify the controlling constituent.
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


def yield_state(cell, step_name, base_strain, strain_direction,
                strain_limit, criterion):
    """Interpolate yield from increments accumulated during one probe."""
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

        # The operational yield criterion uses work increments measured from
        # the base frame of the probe, not values accumulated during preload.
        plastic_increment = np.maximum(
            plastic[:, 1] - plastic[0, 1], 0.0)
        internal_increment = internal[:, 1] - internal[0, 1]
        ratio = (plastic_increment
                 / np.maximum(np.abs(internal_increment), 1.0e-30))

        # First bracket the prescribed dissipation ratio in pseudo-time.
        yielded = np.flatnonzero(ratio >= criterion)
        if yielded.size == 0:
            return None
        upper = int(yielded[0])
        if upper == 0:
            raise RuntimeError(
                "The initial frame satisfies the yield criterion.")
        lower = upper - 1

        # The PEEQ guard prevents an energy-ratio crossing from being accepted
        # before either structural constituent has accumulated plastic strain.
        _, onset_times = constituent_peeq_increment_data(
            cell, step_name, plastic[0, 0], plastic[-1, 0])
        if np.all(np.isnan(onset_times)):
            return None
        plastic_onset = np.nanmin(onset_times)
        if plastic[upper, 0] < plastic_onset - 1.0e-12:
            onset_frame = int(np.searchsorted(
                plastic[:, 0], plastic_onset - 1.0e-12))
            eligible = np.flatnonzero(
                (np.arange(len(ratio)) >= onset_frame)
                & (ratio >= criterion))
            if eligible.size == 0:
                return None
            upper = int(eligible[0])
            lower = upper - 1

        # Interpolate the crossing within the two converged bracket frames.
        denominator = ratio[upper] - ratio[lower]
        if (plastic[upper, 0] <= plastic_onset + 1.0e-12
                and ratio[lower] >= criterion):
            fraction = 1.0
        else:
            fraction = ((criterion - ratio[lower]) / denominator
                        if abs(denominator) > 1.0e-30 else 1.0)
        yield_time = (plastic[lower, 0]
                      + fraction * (plastic[upper, 0] - plastic[lower, 0]))

        # Apply the same interpolation fraction to stress and PEEQ output.
        peeq, _ = constituent_peeq_increment_data(
            cell, step_name, plastic[lower, 0], plastic[upper, 0])
        frame_times = np.asarray(
            [frame.frameValue for frame in step.frames], dtype=float)
        lower_frame = int(np.argmin(
            np.abs(frame_times - plastic[lower, 0])))
        upper_frame = int(np.argmin(
            np.abs(frame_times - plastic[upper, 0])))
        lower_stress = tensor_to_voigt(
            cell.average_stress(step_name, lower_frame), "STRESS")
        upper_stress = tensor_to_voigt(
            cell.average_stress(step_name, upper_frame), "STRESS")
        stress = lower_stress + fraction * (upper_stress - lower_stress)
        strain = base_strain + yield_time * strain_limit * strain_direction
        peeq_yield = peeq[0] + fraction * (peeq[1] - peeq[0])

        # Report which structural constituent first activated the PEEQ guard.
        if np.all(np.isnan(onset_times)):
            controlling_constituent = "NONE"
        elif np.isnan(onset_times[1]) or (
                not np.isnan(onset_times[0])
                and onset_times[0] < onset_times[1] - 1.0e-12):
            controlling_constituent = "SHELL"
        elif np.isnan(onset_times[0]) or (
                onset_times[1] < onset_times[0] - 1.0e-12):
            controlling_constituent = "BEAM"
        else:
            controlling_constituent = "BOTH"

        return {
            "time": yield_time,
            "strain": strain,
            "stress": stress,
            "peeq_shell": peeq_yield[0],
            "peeq_beam": peeq_yield[1],
            "constituent": controlling_constituent,
        }
    finally:
        cell.close_odb_file()


def sample_direction(gmsh_cell, parameters, preload_strain, base_strain,
                     strain_direction, strain_limits, criterion,
                     maximum_strain_increment):
    """Prepare one RVE state and probe a proportional loading direction."""
    preload_strain = np.asarray(preload_strain, dtype=float)
    base_strain = np.asarray(base_strain, dtype=float)
    prepared_state = (np.any(np.abs(preload_strain) > 0.0)
                      or np.any(np.abs(base_strain) > 0.0))

    # Increase the probe range only when the preceding range does not reach
    # yield. Every attempt starts from a fresh RVE with an independent history.
    for strain_limit in strain_limits:
        number_of_increments = int(np.ceil(
            strain_limit / maximum_strain_increment))
        probe_increment = 1.0 / number_of_increments

        def strain_path(time, direction=strain_direction,
                        limit=strain_limit):
            """Return preparation and proportional probing strains."""
            if not prepared_state:
                strain = time * limit * direction
            elif time <= 1.0:
                strain = time * preload_strain
            elif time <= 2.0:
                strain = (preload_strain
                          + (time - 1.0)
                          * (base_strain - preload_strain))
            else:
                strain = (base_strain
                          + (time - 2.0) * limit * direction)
            return voigt_to_tensor(strain, "STRAIN")

        cell = AbaqusRVE(gmsh_cell, custom_types={1: "B31"})
        assign_sections(cell, parameters)

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

        # Applying the complete history at its first step propagates the same
        # periodic macroscopic loading through all subsequent analysis steps.
        cell.apply_periodic_BC(first_step, history)
        cell.model.EmbeddedRegion(
            "EMBEDDED_BEAM_CORE",
            cell.instance.sets["ELSET_BEAM_STRUTS"],
            cell.instance.sets["ELSET_CONTINUUM_MATRIX"])
        cell.run_analysis(
            field_outputs=["PEEQ"],
            history_outputs=["ALLPD", "ALLIE"], numCpus=1)

        # Confirm that returning from the preload to the common base state did
        # not produce additional plastic strain during unloading.
        if prepared_state:
            unload_peeq, _ = constituent_peeq_increment_data(
                cell, unload_step, 0.0, 1.0)
            maximum_unload_peeq = unload_peeq[-1]
            if np.max(maximum_unload_peeq) > 1.0e-10:
                raise RuntimeError(
                    "The common base state was not reached by elastic "
                    "unloading: maximum constituent PEEQ increments are "
                    "{}.".format(maximum_unload_peeq))
            base_stress = tensor_to_voigt(
                cell.average_stress(unload_step), "STRESS")
        else:
            maximum_unload_peeq = np.zeros(2)
            base_stress = np.zeros(6)

        # Extract the first acceptable yield point from the completed probe.
        state = yield_state(
            cell, probe_step, base_strain, strain_direction,
            strain_limit, criterion)
        if state is not None:
            state.update({
                "strain_limit": strain_limit,
                "base_stress": base_stress,
                "unload_peeq_shell": maximum_unload_peeq[0],
                "unload_peeq_beam": maximum_unload_peeq[1],
                # The 3D worker removes this temporary database after the
                # sampled state has been persisted.
                "odb_path": os.path.abspath(cell.job_name + ".odb"),
            })
            return state
    raise RuntimeError(
        "The yield criterion was not reached up to a maximum strain of {}."
        .format(strain_limits[-1]))


# ---------------------------------------------------------------
# Resumable result persistence

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
        except FileExistsError:
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
        with open(path, "r", newline="") as csv_file:
            reader = csv.DictReader(csv_file)
            if reader.fieldnames != ["setting", "value"]:
                raise RuntimeError(
                    "The result configuration '{}' has an unexpected "
                    "header.".format(path))
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
        with os.fdopen(descriptor, "w", newline="") as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(["setting", "value"])
            for key in sorted(expected):
                writer.writerow([key, expected[key]])
        os.replace(temporary_path, path)
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


def write_yield_results(path, records):
    """Write sampled yield points and their common base state."""
    header = (["theta", "yield-time", "strain-limit"]
              + ["epsilon-{}".format(name)
                 for name in ("11", "22", "33", "23", "13", "12")]
              + ["sigma-{}".format(name)
                 for name in ("11", "22", "33", "23", "13", "12")]
              + ["base-sigma-{}".format(name)
                 for name in ("11", "22", "33", "23", "13", "12")]
              + ["max-unload-delta-peeq-shell",
                 "max-unload-delta-peeq-beam",
                 "max-delta-peeq-shell", "max-delta-peeq-beam",
                 "controlling-constituent"])
    with open(path, "w", newline="") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(header)
        for record in records:
            writer.writerow(
                [record["theta"], record["time"],
                 record["strain_limit"]]
                + list(record["strain"])
                + list(record["stress"])
                + list(record["base_stress"])
                + [record["unload_peeq_shell"],
                   record["unload_peeq_beam"],
                   record["peeq_shell"], record["peeq_beam"],
                   record["constituent"]])


def read_yield_results(path):
    """Read completed yield directions for a resumable sampling run."""
    if not os.path.isfile(path):
        return []

    records = []
    with open(path, "r", newline="") as csv_file:
        for row in csv.DictReader(csv_file):
            records.append({
                "theta": float(row["theta"]),
                "time": float(row["yield-time"]),
                "strain_limit": float(row["strain-limit"]),
                "strain": np.asarray([
                    float(row["epsilon-" + name])
                    for name in ("11", "22", "33", "23", "13", "12")]),
                "stress": np.asarray([
                    float(row["sigma-" + name])
                    for name in ("11", "22", "33", "23", "13", "12")]),
                "base_stress": np.asarray([
                    float(row["base-sigma-" + name])
                    for name in ("11", "22", "33", "23", "13", "12")]),
                "unload_peeq_shell": float(
                    row["max-unload-delta-peeq-shell"]),
                "unload_peeq_beam": float(
                    row["max-unload-delta-peeq-beam"]),
                "peeq_shell": float(row["max-delta-peeq-shell"]),
                "peeq_beam": float(row["max-delta-peeq-beam"]),
                "constituent": row["controlling-constituent"],
            })
    return records


# ---------------------------------------------------------------
# Geometry, state definitions, and plane sampling

def build_gmsh_cell(mesh_parameters, model_name="EXAMPLE_3_YIELD"):
    """Build the periodic mesh shared by the yield-surface analyses."""
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
            "result_directory": "results",
        },
        # The evolved surface is probed from the common elastically unloaded
        # state following monotonic epsilon_11 preloading.
        "evolved": {
            "preload_strain": np.asarray(
                [0.03, 0.0, 0.0, 0.0, 0.0, 0.0]),
            "base_strain": np.asarray(
                [0.027, 0.0, 0.0, 0.0, 0.0, 0.0]),
            "strain_limits": (0.003, 0.006, 0.012, 0.024),
            "result_directory": "evolved_results",
        },
    }
    if state_name not in definitions:
        raise ValueError("The yield state must be initial or evolved.")
    return definitions[state_name]


def elastic_calibration(gmsh_cell, parameters, strain_magnitude,
                        result_directory):
    """Return the effective elastic compliance used for stress directions."""
    stiffness_path = os.path.join(
        result_directory, "effective_stiffness.csv")
    compliance_path = os.path.join(
        result_directory, "effective_compliance.csv")

    # Reuse the six elastic calibration analyses across both yield surfaces.
    if os.path.isfile(stiffness_path) and os.path.isfile(compliance_path):
        stiffness = np.loadtxt(stiffness_path, delimiter=",")
        compliance = np.loadtxt(compliance_path, delimiter=",")
        print("Reusing the existing effective elastic calibration.")
    else:
        stiffness = compute_elastic_stiffness(
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
        np.savetxt(stiffness_path, stiffness, delimiter=",")
        np.savetxt(compliance_path, compliance, delimiter=",")
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


def sample_plane(gmsh_cell, parameters, compliance, state_name,
                 preload_strain, base_strain, plane, angles,
                 strain_limits, criterion, maximum_strain_increment,
                 result_directory):
    """Sample one resumable normal-stress-plane section."""
    plane_name, first, second = plane
    result_path = os.path.join(result_directory, plane_name + ".csv")

    # A completed prefix of angular directions can be resumed safely only when
    # it belongs to the same ordered angular grid.
    records = read_yield_results(result_path)
    completed_angles = [record["theta"] for record in records]
    if (len(records) > len(angles)
            or not np.allclose(completed_angles, angles[:len(records)])):
        print("Existing {} results do not match the angular grid; "
              "starting this plane again.".format(plane_name))
        records = []

    print("\nSampling {} {} with {} fixed-increment directions.".format(
        state_name, plane_name, len(angles)))
    if records:
        print("  resuming after {} completed directions".format(
            len(records)))

    for direction_index in range(len(records), len(angles)):
        angle = angles[direction_index]

        # Prescribe a target stress-space direction and use the effective
        # elastic compliance to obtain its proportional strain-space direction.
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
                    strain_direction, strain_limits, criterion,
                    maximum_strain_increment)
                break
            except Exception as error:
                if attempt == 1:
                    raise
                print("    analysis failed once; rebuilding and retrying: "
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

def main():
    """Compute initial and evolved effective yield-surface sections."""

    # Analysis configuration. Set selected_state to "initial", "evolved", or
    # "both". Set selected_plane to one of the names below or leave it as None.
    calibration_only = False
    number_of_directions = 36
    selected_state = "both"
    selected_plane = None
    criterion = 1.0e-4
    maximum_strain_increment = 6.0e-5
    elastic_strain = 1.0e-5
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
    }
    mesh_parameters = {
        "edge_length": 1.0,
        "matrix_elements_per_edge": 8,
        "elements_per_strut": 16,
        "element_order": 1,
        "algorithm_2d": 6,
        "algorithm_3d": 1,
    }

    original_directory = os.getcwd()
    try:
        gmsh_cell = build_gmsh_cell(mesh_parameters)

        # The elastic calibration is shared by the initial and evolved sweeps.
        yield_directory = os.path.join(here, "yield_surface")
        calibration_directory = os.path.join(yield_directory, "results")
        calibration_analysis = os.path.join(
            yield_directory, "calibration_abq")
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
                "elastic_strain": elastic_strain,
                "materials": parameters,
                "mesh": mesh_parameters,
            },
            calibration_results)
        os.chdir(calibration_analysis)
        compliance = elastic_calibration(
            gmsh_cell, parameters, elastic_strain, calibration_directory)
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
                yield_directory, definition["result_directory"])
            analysis_directory = os.path.join(
                yield_directory, state_name + "_abq")
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
                    "analysis": "yield_surface_sections",
                    "criterion": criterion,
                    "directions": number_of_directions,
                    "materials": parameters,
                    "maximum_strain_increment": maximum_strain_increment,
                    "mesh": mesh_parameters,
                    "planes": plane_names,
                    "preload_strain": definition["preload_strain"],
                    "base_strain": definition["base_strain"],
                    "state": state_name,
                    "strain_limits": definition["strain_limits"],
                },
                result_paths)
            os.chdir(analysis_directory)
            np.savetxt(
                os.path.join(result_directory, "base_strain.csv"),
                definition["base_strain"][None, :], delimiter=",")

            print("\n{} yield surface".format(state_name.capitalize()))
            print("  preload strain: {}".format(
                definition["preload_strain"]))
            print("  base strain: {}".format(definition["base_strain"]))
            for plane in stress_planes:
                if (selected_plane is not None
                        and plane[0] != selected_plane):
                    continue
                sample_plane(
                    gmsh_cell, parameters, compliance, state_name,
                    definition["preload_strain"],
                    definition["base_strain"], plane, angles,
                    definition["strain_limits"], criterion,
                    maximum_strain_increment,
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
                base_stresses = np.asarray([
                    record["base_stress"] for record in records])
                if not np.allclose(
                        base_stresses, base_stresses[0],
                        rtol=1.0e-6, atol=1.0e-9):
                    raise RuntimeError(
                        "The independently reproduced base states do not "
                        "have consistent effective stresses.")
                np.savetxt(
                    os.path.join(result_directory, "base_stress.csv"),
                    base_stresses[:1], delimiter=",")

        print("*** Yield-surface sampling finished ***")
    finally:
        os.chdir(original_directory)
        if gmsh.is_initialized():
            gmsh.finalize()


if __name__ == "__main__":
    main()
