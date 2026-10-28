"""Sample sections of the initial effective yield surface.

Every loading direction starts from a fresh, stress-free RVE. Target stress
directions are mapped to proportional strain directions with the numerically
computed effective elastic compliance. Yield is identified when the ratio of
accumulated plastic dissipation to internal energy reaches ``1.0e-4``.

Run with

    abq2024le cae noGUI=yield_surface.py

The default discretization uses 36 independently solved directions over the
complete circle for each stress-plane section. The environment variable
``ABQHOM_YIELD_DIRECTIONS`` can be used for a coarser development run. A
partially completed section is resumed from its CSV file.
"""

import csv
import inspect
import os
import sys
import time

import numpy as np


here = os.path.dirname(os.path.abspath(
    inspect.getframeinfo(inspect.currentframe()).filename))
abqhom_path = os.path.abspath(os.path.join(here, os.pardir,
                                           os.pardir, os.pardir, "src"))
gmsh_path = "C:/Users/jonat/Documents/gmsh/lib"
for path in (here, abqhom_path, gmsh_path):
    if path not in sys.path:
        sys.path.insert(0, path)

import gmsh

from abqhom.gmshmodel import GmshRVE
from abqhom.utils import tensor_to_voigt, voigt_to_tensor
from sandwich import create_geometry

try:
    from abaqusConstants import OFF
    from abqhom.abqmodel import AbaqusRVE, LoadHistory
except ImportError:
    OFF = None
    AbaqusRVE = None
    LoadHistory = None


def assign_materials_and_sections(unit_cell, parameters):
    """Assign the elastic matrix and perfectly plastic structural phases."""
    matrix = unit_cell.model.Material(name="MATRIX")
    matrix.Elastic(table=((parameters["matrix_E"],
                           parameters["matrix_nu"]),))

    shell = unit_cell.model.Material(name="SHELL")
    shell.Elastic(table=((parameters["shell_E"],
                          parameters["shell_nu"]),))
    shell.Plastic(table=((parameters["shell_yield"], 0.0),))

    beam = unit_cell.model.Material(name="BEAM")
    beam.Elastic(table=((parameters["beam_E"],
                         parameters["beam_nu"]),))
    beam.Plastic(table=((parameters["beam_yield"], 0.0),))

    unit_cell.assign_solid_section("CONTINUUM_MATRIX", "MATRIX")
    unit_cell.assign_shell_section(
        "SHELL_FACES", "SHELL", parameters["shell_thickness"],
        num_integration_points=5, pre_integrate=False)
    unit_cell.assign_beam_section(
        "BEAM_STRUTS", "BEAM", parameters["beam_radius"],
        pre_integrate=False)


def create_unit_cell(gmsh_unit_cell, parameters):
    """Return a fresh Abaqus model with materials and sections assigned."""
    unit_cell = AbaqusRVE(
        gmsh_unit_cell, custom_types={1: "B31"}, remove_secondary=True)
    assign_materials_and_sections(unit_cell, parameters)
    return unit_cell


def add_embedded_core(unit_cell):
    """Embed the complete beam core after periodic constraints are defined."""
    unit_cell.model.EmbeddedRegion(
        "EMBEDDED_BEAM_CORE",
        unit_cell.instance.sets["ELSET_BEAM_STRUTS"],
        unit_cell.instance.sets["ELSET_CONTINUUM_MATRIX"])


def elastic_stiffness(gmsh_unit_cell, parameters, strain_magnitude):
    """Compute the complete effective elastic stiffness in Voigt notation."""
    stiffness = np.zeros((6, 6))
    for component in range(6):
        print("Elastic stiffness case {}/6".format(component + 1))
        sys.stdout.flush()

        def strain_path(time, index=component):
            strain = np.zeros(6)
            strain[index] = time*strain_magnitude
            return voigt_to_tensor(strain, "STRAIN")

        history = LoadHistory(np.asarray([0.0, 1.0]), strain_path)
        unit_cell = create_unit_cell(gmsh_unit_cell, parameters)
        step_name = "ELASTIC_LOAD"
        unit_cell.model.StaticStep(
            name=step_name, previous="Initial", nlgeom=OFF,
            timePeriod=1.0, maxNumInc=1,
            initialInc=1.0, minInc=1.0, maxInc=1.0)
        unit_cell.apply_periodic_BC(step_name, history)
        add_embedded_core(unit_cell)
        run_analysis_and_wait(unit_cell, numCpus=1)
        stiffness[:, component] = (
            tensor_to_voigt(unit_cell.average_stress(step_name), "STRESS")
            / strain_magnitude)

    return stiffness


def phase_peeq_data(unit_cell, step_name, lower_time, upper_time,
                    tolerance=1.0e-12):
    """Return phase-wise PEEQ values and first-yield times.

    Because PEEQ is non-decreasing, the first frame with nonzero PEEQ is
    located by a binary search rather than by reading every output frame.
    """
    close_odb = unit_cell.odb is None
    unit_cell.open_odb_file()
    try:
        step = unit_cell.odb.steps[step_name]
        instance = unit_cell.odb.rootAssembly.instances[unit_cell.model_name]
        phase_sets = (
            instance.elementSets["ELSET_SHELL_FACES"],
            instance.elementSets["ELSET_BEAM_STRUTS"],
        )
        frames = step.frames
        frame_times = np.asarray(
            [frame.frameValue for frame in frames], dtype=float)
        lower_frame = int(np.argmin(np.abs(frame_times - lower_time)))
        upper_frame = int(np.argmin(np.abs(frame_times - upper_time)))
        cache = {}

        def phase_maximum(frame_index, phase):
            key = (frame_index, phase)
            if key not in cache:
                peeq = frames[frame_index].fieldOutputs["PEEQ"]
                values = peeq.getSubset(region=phase_sets[phase]).values
                cache[key] = max(float(value.data) for value in values)
            return cache[key]

        onset_times = [np.nan, np.nan]
        for phase in range(len(phase_sets)):
            if phase_maximum(len(frames) - 1, phase) <= tolerance:
                continue
            left = 0
            right = len(frames) - 1
            while left < right:
                middle = (left + right)//2
                if phase_maximum(middle, phase) > tolerance:
                    right = middle
                else:
                    left = middle + 1
            onset_times[phase] = frame_times[left]

        peeq = np.asarray([
            [phase_maximum(lower_frame, phase)
             for phase in range(len(phase_sets))],
            [phase_maximum(upper_frame, phase)
             for phase in range(len(phase_sets))],
        ])
    finally:
        if close_odb:
            unit_cell.close_odb_file()
    return peeq, np.asarray(onset_times)


def yield_state(unit_cell, step_name, strain_direction, strain_limit,
                criterion):
    """Interpolate the effective state at the selected dissipation ratio."""
    unit_cell.open_odb_file()
    try:
        step = unit_cell.odb.steps[step_name]
        history_region = step.historyRegions["Assembly ASSEMBLY"]
        plastic = np.asarray(
            history_region.historyOutputs["ALLPD"].data, dtype=float)
        internal = np.asarray(
            history_region.historyOutputs["ALLIE"].data, dtype=float)
        if not np.allclose(plastic[:, 0], internal[:, 0]):
            raise RuntimeError(
                "ALLPD and ALLIE are not available at identical times.")
        plastic_increment = np.maximum(
            plastic[:, 1] - plastic[0, 1], 0.0)
        internal_increment = internal[:, 1] - internal[0, 1]
        ratio = (plastic_increment
                 / np.maximum(np.abs(internal_increment), 1.0e-30))

        yielded = np.flatnonzero(ratio >= criterion)
        if yielded.size == 0:
            return None
        upper = int(yielded[0])
        if upper == 0:
            raise RuntimeError(
                "The initial frame satisfies the yield criterion.")
        lower = upper - 1
        peeq, onset_times = phase_peeq_data(
            unit_cell, step_name, plastic[lower, 0], plastic[upper, 0])
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

        denominator = ratio[upper] - ratio[lower]
        if (plastic[upper, 0] <= plastic_onset + 1.0e-12
                and ratio[lower] >= criterion):
            fraction = 1.0
        else:
            fraction = ((criterion - ratio[lower])/denominator
                        if abs(denominator) > 1.0e-30 else 1.0)
        yield_time = (plastic[lower, 0]
                      + fraction*(plastic[upper, 0] - plastic[lower, 0]))

        peeq, onset_times = phase_peeq_data(
            unit_cell, step_name, plastic[lower, 0], plastic[upper, 0])
        frame_times = np.asarray(
            [frame.frameValue for frame in step.frames], dtype=float)
        lower_frame = int(np.argmin(
            np.abs(frame_times - plastic[lower, 0])))
        upper_frame = int(np.argmin(
            np.abs(frame_times - plastic[upper, 0])))
        lower_stress = tensor_to_voigt(
            unit_cell.average_stress(step_name, lower_frame), "STRESS")
        upper_stress = tensor_to_voigt(
            unit_cell.average_stress(step_name, upper_frame), "STRESS")
        stress = lower_stress + fraction*(upper_stress - lower_stress)
        strain = yield_time*strain_limit*strain_direction
        peeq_yield = peeq[0] + fraction*(peeq[1] - peeq[0])

        if np.all(np.isnan(onset_times)):
            controlling_phase = "NONE"
        elif np.isnan(onset_times[1]) or (
                not np.isnan(onset_times[0])
                and onset_times[0] < onset_times[1] - 1.0e-12):
            controlling_phase = "SHELL"
        elif np.isnan(onset_times[0]) or (
                onset_times[1] < onset_times[0] - 1.0e-12):
            controlling_phase = "BEAM"
        else:
            controlling_phase = "BOTH"

        return {
            "time": yield_time,
            "strain": strain,
            "stress": stress,
            "peeq_shell": peeq_yield[0],
            "peeq_beam": peeq_yield[1],
            "phase": controlling_phase,
        }
    finally:
        unit_cell.close_odb_file()


def sample_direction(gmsh_unit_cell, parameters, strain_direction,
                     strain_limits, criterion, increment_size):
    """Run one virgin RVE and return the first detected yield state."""
    number_of_increments = int(round(1.0/increment_size))
    if not np.isclose(number_of_increments*increment_size, 1.0):
        raise ValueError(
            "The fixed increment size must divide the unit step exactly.")

    for strain_limit in strain_limits:
        def strain_path(time, direction=strain_direction,
                        limit=strain_limit):
            return voigt_to_tensor(time*limit*direction, "STRAIN")

        history = LoadHistory(np.asarray([0.0, 1.0]), strain_path)
        unit_cell = create_unit_cell(gmsh_unit_cell, parameters)
        step_name = "PROPORTIONAL_LOAD"
        unit_cell.model.StaticStep(
            name=step_name, previous="Initial", nlgeom=OFF,
            timePeriod=1.0, maxNumInc=number_of_increments,
            initialInc=increment_size, minInc=increment_size,
            maxInc=increment_size)
        unit_cell.apply_periodic_BC(step_name, history)
        add_embedded_core(unit_cell)
        run_analysis_and_wait(
            unit_cell,
            field_outputs=["PEEQ"],
            history_outputs=["ALLPD", "ALLIE"], numCpus=1)
        state = yield_state(
            unit_cell, step_name, strain_direction, strain_limit,
            criterion)
        if state is not None:
            state["strain_limit"] = strain_limit
            return state
    raise RuntimeError(
        "The yield criterion was not reached up to a maximum strain of {}."
        .format(strain_limits[-1]))


def run_analysis_and_wait(unit_cell, timeout=3600.0, **kwargs):
    """Run an analysis and wait until its new ODB is ready for reading."""
    odb_path = os.path.abspath(
        "JOB_{}.odb".format(unit_cell.model_name))
    lock_path = os.path.splitext(odb_path)[0] + ".lck"
    log_path = os.path.splitext(odb_path)[0] + ".log"
    unit_cell.close_odb_file()
    if os.path.isfile(lock_path):
        raise RuntimeError(
            "The working directory contains the active lock {}."
            .format(lock_path))
    if os.path.isfile(odb_path):
        os.remove(odb_path)

    unit_cell.run_analysis(**kwargs)
    start_time = time.time()
    while not os.path.isfile(odb_path) or os.path.isfile(lock_path):
        if time.time() - start_time > timeout:
            raise RuntimeError(
                "Timed out while waiting for {}.".format(odb_path))
        time.sleep(1.0)
    if os.path.isfile(log_path):
        with open(log_path, "r") as log_file:
            log_text = log_file.read()
        if "Abaqus JOB {} COMPLETED".format(
                "JOB_{}".format(unit_cell.model_name)) not in log_text:
            raise RuntimeError(
                "Abaqus analysis did not complete successfully; see {}."
                .format(log_path))
    return odb_path


def write_yield_results(path, records):
    """Write sampled yield points and their complete effective states."""
    header = (["theta", "yield_time", "strain-limit"]
              + ["epsilon-{}".format(name)
                 for name in ("11", "22", "33", "23", "13", "12")]
              + ["sigma-{}".format(name)
                 for name in ("11", "22", "33", "23", "13", "12")]
              + ["max-peeq-shell", "max-peeq-beam",
                 "controlling-phase"])
    with open(path, "w", newline="") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(header)
        for record in records:
            writer.writerow(
                [record["theta"], record["time"],
                 record["strain_limit"]]
                + list(record["strain"])
                + list(record["stress"])
                + [record["peeq_shell"], record["peeq_beam"],
                   record["phase"]])


def read_yield_results(path):
    """Read completed yield directions for a resumable sampling run."""
    if not os.path.isfile(path):
        return []

    records = []
    with open(path, "r", newline="") as csv_file:
        rows = csv.reader(csv_file)
        next(rows)
        for row in rows:
            records.append({
                "theta": float(row[0]),
                "time": float(row[1]),
                "strain_limit": float(row[2]),
                "strain": np.asarray([float(value) for value in row[3:9]]),
                "stress": np.asarray([float(value) for value in row[9:15]]),
                "peeq_shell": float(row[15]),
                "peeq_beam": float(row[16]),
                "phase": row[17],
            })
    return records


def build_gmsh_unit_cell(model_name="EXAMPLE_3_YIELD"):
    """Build the periodic mesh shared by the yield-surface analyses."""
    edge_length = 1.0
    create_geometry(model_name, edge_length, elements_per_strut=16)
    gmsh_unit_cell = GmshRVE(model_name)
    gmsh_unit_cell.add_periodic_mesh_constraints()
    gmsh_unit_cell.generate_mesh(
        element_order=1, algorithm_2d=6, algorithm_3d=1,
        element_size=edge_length/8.0)
    gmsh_unit_cell.build_node_sets()
    gmsh_unit_cell.build_element_sets()
    return gmsh_unit_cell


def constituent_parameters():
    """Return the parameters used by the cyclic and yield-surface runs."""
    return {
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


def main():
    """Compute two sections of the initial effective yield surface."""
    gmsh_unit_cell = build_gmsh_unit_cell()
    if AbaqusRVE is None:
        gmsh.finalize()
        print("Abaqus modules are unavailable; mesh generation is complete.")
        return

    parameters = constituent_parameters()
    yield_directory = os.path.join(here, "yield_surface")
    analysis_name = os.environ.get(
        "ABQHOM_YIELD_ANALYSIS_DIRECTORY", "initial_abq")
    analysis_directory = os.path.join(yield_directory, analysis_name)
    result_directory = os.path.join(yield_directory, "results")
    for directory in (analysis_directory, result_directory):
        if not os.path.isdir(directory):
            os.makedirs(directory)
    os.chdir(analysis_directory)

    stiffness_path = os.path.join(
        result_directory, "effective_stiffness.csv")
    compliance_path = os.path.join(
        result_directory, "effective_compliance.csv")
    if os.path.isfile(stiffness_path) and os.path.isfile(compliance_path):
        stiffness = np.loadtxt(stiffness_path, delimiter=",")
        compliance = np.loadtxt(compliance_path, delimiter=",")
        print("Reusing the existing effective elastic calibration.")
    else:
        elastic_strain = 1.0e-5
        stiffness = elastic_stiffness(
            gmsh_unit_cell, parameters, elastic_strain)
        compliance = None
    stiffness_symmetric = 0.5*(stiffness + stiffness.T)
    symmetry_residual = (np.linalg.norm(stiffness - stiffness.T)
                         / np.linalg.norm(stiffness))
    eigenvalues = np.linalg.eigvalsh(stiffness_symmetric)
    if np.min(eigenvalues) <= 0.0:
        raise RuntimeError("The effective elastic stiffness is not positive.")
    if compliance is None:
        compliance = np.linalg.inv(stiffness_symmetric)
        np.savetxt(stiffness_path, stiffness, delimiter=",")
        np.savetxt(compliance_path, compliance, delimiter=",")
    print("Elastic stiffness symmetry residual: {:.3e}".format(
        symmetry_residual))
    print("Elastic stiffness eigenvalues: {}".format(eigenvalues))
    sys.stdout.flush()
    if os.environ.get("ABQHOM_YIELD_CALIBRATION_ONLY", "0") == "1":
        gmsh.finalize()
        print("*** Elastic calibration finished ***")
        return

    number_of_directions = int(os.environ.get(
        "ABQHOM_YIELD_DIRECTIONS", "36"))
    if number_of_directions < 4:
        raise ValueError(
            "At least four yield-surface directions are required.")
    angles = np.linspace(
        0.0, 2.0*np.pi, number_of_directions, endpoint=False)
    yield_criterion = 1.0e-4
    increment_size = 0.002
    strain_limits = (0.03, 0.06)
    stress_planes = (
        ("sigma11_sigma22", 0, 1),
        ("sigma11_sigma33", 0, 2),
    )
    selected_plane = os.environ.get("ABQHOM_YIELD_PLANE")

    for plane_name, first, second in stress_planes:
        if selected_plane is not None and plane_name != selected_plane:
            continue
        result_path = os.path.join(result_directory, plane_name + ".csv")
        records = read_yield_results(result_path)
        completed_angles = [record["theta"] for record in records]
        if (len(records) > number_of_directions
                or not np.allclose(
                    completed_angles, angles[:len(records)])):
            print("Existing {} results do not match the angular grid; "
                  "starting this plane again.".format(plane_name))
            records = []
        print("\nSampling {} with {} directions and fixed increment "
              "Delta t = {}.".format(
                  plane_name, number_of_directions, increment_size))
        if records:
            print("  resuming after {} completed directions".format(
                len(records)))
        direction_indices = range(len(records), number_of_directions)
        for direction_index in direction_indices:
            angle = angles[direction_index]
            stress_direction = np.zeros(6)
            stress_direction[first] = np.cos(angle)
            stress_direction[second] = np.sin(angle)
            strain_direction = compliance.dot(stress_direction)
            strain_direction /= np.max(np.abs(strain_direction))

            print("  direction {:02d}/{:02d}: theta = {:7.2f} deg".format(
                direction_index + 1, number_of_directions,
                np.degrees(angle)))
            sys.stdout.flush()
            for attempt in range(2):
                try:
                    state = sample_direction(
                        gmsh_unit_cell, parameters, strain_direction,
                        strain_limits, yield_criterion, increment_size)
                    break
                except Exception as error:
                    if attempt == 1:
                        raise
                    print("    analysis failed once; rebuilding and retrying: "
                          "{}".format(error))
                    sys.stdout.flush()
            state["theta"] = angle
            if direction_index < len(records):
                records[direction_index] = state
            else:
                records.append(state)
            write_yield_results(result_path, records)

        write_yield_results(result_path, records)

        print("Completed {}.".format(plane_name))
        sys.stdout.flush()

    gmsh.finalize()
    print("*** Yield-surface sampling finished ***")


if __name__ == "__main__":
    main()

