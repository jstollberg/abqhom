"""Sample the initial or evolved yield surface in normal-stress space.

The sampling directions are the 642 vertices of a level-three icosphere. Both
the virgin and evolved surfaces are sampled independently on the complete
sphere.

Use ``ABQHOM_YIELD_3D_STATE`` to select ``initial`` or ``final``. Independent
workers are selected with ``ABQHOM_YIELD_3D_WORKER`` and
``ABQHOM_YIELD_3D_WORKERS``. Each worker writes an independent, resumable CSV
file. Set ``ABQHOM_YIELD_3D_MERGE=1`` after all workers finish to assemble the
complete surface and its triangular connectivity.
"""

import csv
import inspect
import os
import sys

import numpy as np


here = os.path.dirname(os.path.abspath(
    inspect.getframeinfo(inspect.currentframe()).filename))
if here not in sys.path:
    sys.path.insert(0, here)

from yield_surface import (AbaqusRVE, LoadHistory, OFF, add_embedded_core,
                           build_gmsh_unit_cell, constituent_parameters,
                           create_unit_cell, gmsh, run_analysis_and_wait,
                           tensor_to_voigt,
                           voigt_to_tensor, yield_state)
from final_yield_surface import (base_strain, complete_strain_path,
                                 final_yield_state,
                                 phase_peeq_increment_data)

component_names = ("11", "22", "33", "23", "13", "12")
header = (["direction-index", "n1", "n2", "n3", "yield-time",
           "strain-limit"]
          + ["epsilon-{}".format(name) for name in component_names]
          + ["sigma-{}".format(name) for name in component_names]
          + ["base-sigma-{}".format(name) for name in component_names]
          + ["max-unload-delta-peeq-shell",
             "max-unload-delta-peeq-beam",
             "max-delta-peeq-shell", "max-delta-peeq-beam",
             "controlling-phase"])


def icosphere(subdivisions=3):
    """Return unit vertices and triangular faces of an icosphere."""
    golden_ratio = 0.5*(1.0 + np.sqrt(5.0))
    vertices = [
        (-1.0, golden_ratio, 0.0), (1.0, golden_ratio, 0.0),
        (-1.0, -golden_ratio, 0.0), (1.0, -golden_ratio, 0.0),
        (0.0, -1.0, golden_ratio), (0.0, 1.0, golden_ratio),
        (0.0, -1.0, -golden_ratio), (0.0, 1.0, -golden_ratio),
        (golden_ratio, 0.0, -1.0), (golden_ratio, 0.0, 1.0),
        (-golden_ratio, 0.0, -1.0), (-golden_ratio, 0.0, 1.0),
    ]
    faces = [
        (0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11),
        (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
        (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9),
        (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1),
    ]
    vertices = [np.asarray(vertex, dtype=float)
                / np.linalg.norm(vertex) for vertex in vertices]

    for _ in range(subdivisions):
        midpoint_cache = {}

        def midpoint(first, second):
            edge = tuple(sorted((first, second)))
            if edge not in midpoint_cache:
                vertex = vertices[first] + vertices[second]
                vertex /= np.linalg.norm(vertex)
                midpoint_cache[edge] = len(vertices)
                vertices.append(vertex)
            return midpoint_cache[edge]

        refined_faces = []
        for first, second, third in faces:
            first_second = midpoint(first, second)
            second_third = midpoint(second, third)
            third_first = midpoint(third, first)
            refined_faces.extend((
                (first, first_second, third_first),
                (second, second_third, first_second),
                (third, third_first, second_third),
                (first_second, second_third, third_first),
            ))
        faces = refined_faces

    return np.asarray(vertices), np.asarray(faces, dtype=int)


def remove_odb(path):
    """Remove a closed temporary ODB after its result has been persisted."""
    if path is not None and os.path.isfile(path):
        os.remove(path)


def sample_initial_direction(gmsh_unit_cell, parameters, strain_direction,
                             strain_limits, criterion,
                             maximum_strain_increment):
    """Return a virgin-state yield point and the temporary ODB path."""
    for strain_limit in strain_limits:
        number_of_increments = int(np.ceil(
            strain_limit/maximum_strain_increment))
        increment = 1.0/number_of_increments

        def strain_path(time, direction=strain_direction,
                        limit=strain_limit):
            return voigt_to_tensor(time*limit*direction, "STRAIN")

        unit_cell = create_unit_cell(gmsh_unit_cell, parameters)
        step_name = "PROPORTIONAL_LOAD"
        unit_cell.model.StaticStep(
            name=step_name, previous="Initial", nlgeom=OFF,
            timePeriod=1.0, maxNumInc=number_of_increments,
            initialInc=increment, minInc=increment, maxInc=increment)
        unit_cell.apply_periodic_BC(
            step_name, LoadHistory(np.asarray([0.0, 1.0]), strain_path))
        add_embedded_core(unit_cell)
        odb_path = run_analysis_and_wait(
            unit_cell,
            field_outputs=["PEEQ"],
            history_outputs=["ALLPD", "ALLIE"], numCpus=1)
        state = yield_state(
            unit_cell, step_name, strain_direction, strain_limit,
            criterion)
        if state is not None:
            state.update({
                "strain_limit": strain_limit,
                "base_stress": np.zeros(6),
                "unload_peeq_shell": 0.0,
                "unload_peeq_beam": 0.0,
            })
            return state, odb_path
        remove_odb(odb_path)
    raise RuntimeError(
        "The initial yield criterion was not reached up to strain {}."
        .format(strain_limits[-1]))


def sample_final_direction(gmsh_unit_cell, parameters, strain_direction,
                           strain_limits, criterion,
                           maximum_strain_increment):
    """Return an evolved-state yield point and the temporary ODB path."""
    for strain_limit in strain_limits:
        number_of_increments = int(np.ceil(
            strain_limit/maximum_strain_increment))
        increment = 1.0/number_of_increments

        def strain_path(time, direction=strain_direction,
                        limit=strain_limit):
            return complete_strain_path(time, direction, limit)

        unit_cell = create_unit_cell(gmsh_unit_cell, parameters)
        preload_step = "MONOTONIC_PRELOAD"
        unload_step = "ELASTIC_UNLOAD"
        probe_step = "PROBE"
        unit_cell.model.StaticStep(
            name=preload_step, previous="Initial", nlgeom=OFF,
            timePeriod=1.0, maxNumInc=100,
            initialInc=0.01, minInc=0.01, maxInc=0.01)
        unit_cell.model.StaticStep(
            name=unload_step, previous=preload_step, nlgeom=OFF,
            timePeriod=1.0, maxNumInc=10,
            initialInc=0.1, minInc=0.1, maxInc=0.1)
        unit_cell.model.StaticStep(
            name=probe_step, previous=unload_step, nlgeom=OFF,
            timePeriod=1.0, maxNumInc=number_of_increments,
            initialInc=increment, minInc=increment, maxInc=increment)
        history = LoadHistory(
            np.asarray([0.0, 1.0, 2.0, 3.0]), strain_path)
        unit_cell.apply_periodic_BC(preload_step, history)
        add_embedded_core(unit_cell)
        odb_path = run_analysis_and_wait(
            unit_cell,
            field_outputs=["PEEQ"],
            history_outputs=["ALLPD", "ALLIE"], numCpus=1)

        unload_peeq, _ = phase_peeq_increment_data(
            unit_cell, unload_step, 0.0, 1.0)
        maximum_unload_peeq = unload_peeq[-1]
        if np.max(maximum_unload_peeq) > 1.0e-10:
            remove_odb(odb_path)
            raise RuntimeError(
                "The common probing state was not reached elastically: {}."
                .format(maximum_unload_peeq))
        base_stress = tensor_to_voigt(
            unit_cell.average_stress(unload_step), "STRESS")
        state = final_yield_state(
            unit_cell, probe_step, base_strain, strain_direction,
            strain_limit, criterion)
        if state is not None:
            state.update({
                "strain_limit": strain_limit,
                "base_stress": base_stress,
                "unload_peeq_shell": maximum_unload_peeq[0],
                "unload_peeq_beam": maximum_unload_peeq[1],
            })
            return state, odb_path
        remove_odb(odb_path)
    raise RuntimeError(
        "The evolved yield criterion was not reached up to strain {}."
        .format(strain_limits[-1]))


def record_row(direction_index, direction, state):
    """Return one CSV row for a sampled state."""
    return ([direction_index] + list(direction)
            + [state["time"], state["strain_limit"]]
            + list(state["strain"]) + list(state["stress"])
            + list(state["base_stress"])
            + [state["unload_peeq_shell"],
               state["unload_peeq_beam"],
               state["peeq_shell"], state["peeq_beam"],
               state["phase"]])


def read_part(path):
    """Return records from one worker file keyed by direction index."""
    if not os.path.isfile(path):
        return {}
    with open(path, newline="") as csv_file:
        rows = list(csv.DictReader(csv_file))
    return {int(row["direction-index"]): row for row in rows}


def append_row(path, row):
    """Append one completed result to a worker CSV."""
    exists = os.path.isfile(path)
    with open(path, "a", newline="") as csv_file:
        writer = csv.writer(csv_file)
        if not exists:
            writer.writerow(header)
        writer.writerow(row)


def write_complete_surface(state_name, root_directory,
                           directions, faces, workers):
    """Merge worker files and write the complete ordered surface."""
    part_directory = os.path.join(root_directory, "parts")
    records = {}
    for worker in range(workers):
        path = os.path.join(
            part_directory, "worker_{:02d}.csv".format(worker))
        worker_records = read_part(path)
        invalid = [index for index in worker_records
                   if index % workers != worker]
        if invalid:
            raise RuntimeError(
                "Worker {} contains directions assigned to another worker."
                .format(worker))
        overlap = set(records).intersection(worker_records)
        if overlap:
            raise RuntimeError(
                "Directions occur in more than one worker file: {}."
                .format(sorted(overlap)))
        records.update(worker_records)

    missing = sorted(set(range(len(directions))).difference(records))
    if missing:
        raise RuntimeError(
            "The {} sampling is missing {} directions."
            .format(state_name, len(missing)))

    output = os.path.join(root_directory, "yield_surface.csv")
    with open(output, "w", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=header)
        writer.writeheader()
        for index in range(len(directions)):
            writer.writerow(records[index])

    np.savetxt(os.path.join(root_directory, "directions.csv"),
               np.column_stack((np.arange(len(directions)), directions)),
               delimiter=",", header="direction-index,n1,n2,n3",
               comments="", fmt=["%d", "%.17e", "%.17e", "%.17e"])
    np.savetxt(os.path.join(root_directory, "faces.csv"), faces,
               delimiter=",", header="vertex-1,vertex-2,vertex-3",
               comments="", fmt="%d")
    print("Wrote {} complete directions and {} faces to {}.".format(
        len(directions), len(faces), output))


def main():
    """Run one sampling worker or merge all completed worker files."""
    state_name = os.environ.get(
        "ABQHOM_YIELD_3D_STATE", "initial").lower()
    if state_name not in ("initial", "final"):
        raise ValueError("The 3D yield state must be initial or final.")
    worker = int(os.environ.get("ABQHOM_YIELD_3D_WORKER", "0"))
    workers = int(os.environ.get("ABQHOM_YIELD_3D_WORKERS", "1"))
    if worker < 0 or worker >= workers:
        raise ValueError("The worker index must lie in [0, workers).")

    directions, faces = icosphere(subdivisions=3)
    root_directory = os.path.join(
        here, "yield_surface", "surface_3d", state_name)
    part_directory = os.path.join(root_directory, "parts")
    try:
        os.makedirs(part_directory)
    except OSError:
        if not os.path.isdir(part_directory):
            raise

    if os.environ.get("ABQHOM_YIELD_3D_MERGE") == "1":
        write_complete_surface(
            state_name, root_directory, directions, faces,
            workers)
        return

    if AbaqusRVE is None:
        raise RuntimeError("The sampling must be run with Abaqus Python.")
    analysis_directory = os.path.join(
        root_directory, "abq_worker_{:02d}".format(worker))
    try:
        os.makedirs(analysis_directory)
    except OSError:
        if not os.path.isdir(analysis_directory):
            raise
    os.chdir(analysis_directory)

    part_path = os.path.join(
        part_directory, "worker_{:02d}.csv".format(worker))
    completed = read_part(part_path)
    if state_name == "initial":
        strain_limits = (0.006, 0.012, 0.024)
    else:
        strain_limits = (0.003, 0.006, 0.012, 0.024)
    selected = list(range(worker, len(directions), workers))

    gmsh_unit_cell = build_gmsh_unit_cell(
        "EXAMPLE_3_YIELD3D_{}_{:02d}".format(
            state_name.upper(), worker))
    parameters = constituent_parameters()
    compliance = np.loadtxt(os.path.join(
        here, "yield_surface", "results", "effective_compliance.csv"),
        delimiter=",")
    criterion = 1.0e-4
    maximum_strain_increment = 6.0e-5
    sampler = (sample_initial_direction if state_name == "initial"
               else sample_final_direction)

    print("{} worker {}/{}: {} assigned directions, {} complete."
          .format(state_name, worker + 1, workers, len(selected),
                  len(set(selected).intersection(completed))))
    sys.stdout.flush()
    for position, direction_index in enumerate(selected, start=1):
        if direction_index in completed:
            continue
        direction = directions[direction_index]
        stress_direction = np.zeros(6)
        stress_direction[:3] = direction
        strain_direction = compliance.dot(stress_direction)
        strain_direction /= np.max(np.abs(strain_direction))
        print("  worker {}: {}/{} direction {} ({:+.4f}, {:+.4f}, "
              "{:+.4f})".format(
                  worker, position, len(selected), direction_index,
                  direction[0], direction[1], direction[2]))
        sys.stdout.flush()

        for attempt in range(2):
            try:
                state, odb_path = sampler(
                    gmsh_unit_cell, parameters, strain_direction,
                    strain_limits, criterion, maximum_strain_increment)
                break
            except Exception as error:
                if attempt == 1:
                    raise
                print("    analysis failed once; retrying: {}".format(error))
                sys.stdout.flush()
        append_row(part_path, record_row(direction_index, direction, state))
        remove_odb(odb_path)

    gmsh.finalize()
    print("*** {} worker {} finished ***".format(state_name, worker))


if __name__ == "__main__":
    main()
