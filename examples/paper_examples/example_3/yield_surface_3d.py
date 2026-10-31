"""Sample the initial or evolved yield surface in normal-stress space.

The sampling directions are the 642 vertices of a level-three icosphere. Both
the initial and evolved surfaces are sampled independently on the complete
sphere.

Select the surface, worker index, worker count, and merge operation in the
configuration block of ``main``. Each worker writes an independent, resumable
CSV file. After all workers finish, the merge operation assembles the complete
surface and its triangular connectivity.
"""

import csv
import inspect
import os
import sys

import gmsh
import numpy as np

# ---------------------------------------------------------------
# Locate the neighboring planar-section script. Under
# 'abaqus cae noGUI=yield_surface_3d.py' __file__ is undefined, so the
# script location is obtained via inspect instead.
here = os.path.dirname(os.path.abspath(
    inspect.getframeinfo(inspect.currentframe()).filename))
if here not in sys.path:
    sys.path.insert(0, here)

from yield_surface_sections import (build_gmsh_cell, elastic_calibration,
                                    sample_direction, state_definition,
                                    validate_result_configuration)

COMPONENT_NAMES = ("11", "22", "33", "23", "13", "12")
CSV_HEADER = (["direction-index", "n1", "n2", "n3", "yield-time",
               "strain-limit"]
              + ["epsilon-{}".format(name) for name in COMPONENT_NAMES]
              + ["sigma-{}".format(name) for name in COMPONENT_NAMES]
              + ["base-sigma-{}".format(name) for name in COMPONENT_NAMES]
              + ["max-unload-delta-peeq-shell",
                 "max-unload-delta-peeq-beam",
                 "max-delta-peeq-shell", "max-delta-peeq-beam",
                 "controlling-constituent"])


# ---------------------------------------------------------------
# Direction grid

def icosphere(subdivisions=3):
    """Return unit vertices and triangular faces of an icosphere."""
    if subdivisions < 0:
        raise ValueError("The number of icosphere subdivisions is negative.")

    # Start from the twelve normalized vertices and twenty faces of a regular
    # icosahedron. Every refinement splits one triangle into four triangles.
    golden_ratio = 0.5 * (1.0 + np.sqrt(5.0))
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
            """Return the cached midpoint vertex of an icosphere edge."""
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


# ---------------------------------------------------------------
# Worker persistence and result assembly

def remove_output_database(path):
    """Remove a closed temporary ODB after its result has been persisted."""
    if path is not None and os.path.isfile(path):
        os.remove(path)


def result_row(direction_index, direction, state):
    """Return one CSV row for a sampled state."""
    return ([direction_index] + list(direction)
            + [state["time"], state["strain_limit"]]
            + list(state["strain"]) + list(state["stress"])
            + list(state["base_stress"])
            + [state["unload_peeq_shell"],
               state["unload_peeq_beam"],
               state["peeq_shell"], state["peeq_beam"],
               state["constituent"]])


def read_worker_results(path):
    """Return records from one worker file keyed by direction index."""
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return {}

    with open(path, newline="") as csv_file:
        reader = csv.DictReader(csv_file)
        if reader.fieldnames != CSV_HEADER:
            raise RuntimeError(
                "Worker result file '{}' has an unexpected header."
                .format(path))
        rows = list(reader)

    records = {}
    for row in rows:
        index = int(row["direction-index"])
        if index in records:
            raise RuntimeError(
                "Direction {} occurs more than once in '{}'."
                .format(index, path))
        records[index] = row
    return records


def validate_worker_directions(records, directions):
    """Validate stored direction vectors against the current sampling grid."""
    for index, record in records.items():
        direction = np.asarray([
            float(record["n1"]), float(record["n2"]), float(record["n3"])
        ])
        if not np.allclose(
                direction, directions[index], rtol=1.0e-12, atol=1.0e-14):
            raise RuntimeError(
                "Stored direction {} does not match the current icosphere."
                .format(index))


def append_worker_result(path, row):
    """Append one completed result to a worker CSV."""
    if len(row) != len(CSV_HEADER):
        raise ValueError("The worker result does not match the CSV header.")

    has_content = os.path.isfile(path) and os.path.getsize(path) > 0
    with open(path, "a", newline="") as csv_file:
        writer = csv.writer(csv_file)
        if not has_content:
            writer.writerow(CSV_HEADER)
        writer.writerow(row)


def merge_worker_results(state_name, root_directory, directions, faces,
                         number_of_workers):
    """Merge worker files and write the complete ordered surface."""
    part_directory = os.path.join(root_directory, "parts")
    records = {}
    valid_indices = set(range(len(directions)))

    # Modulo assignment gives every worker a deterministic, disjoint subset.
    for worker_index in range(number_of_workers):
        path = os.path.join(
            part_directory, "worker_{:02d}.csv".format(worker_index))
        worker_records = read_worker_results(path)
        invalid = [
            index for index in worker_records
            if (index not in valid_indices
                or index % number_of_workers != worker_index)]
        if invalid:
            raise RuntimeError(
                "Worker {} contains invalid or incorrectly assigned "
                "directions: {}.".format(worker_index, sorted(invalid)))
        validate_worker_directions(worker_records, directions)
        overlap = set(records).intersection(worker_records)
        if overlap:
            raise RuntimeError(
                "Directions occur in more than one worker file: {}."
                .format(sorted(overlap)))
        records.update(worker_records)

    missing = sorted(valid_indices.difference(records))
    if missing:
        raise RuntimeError(
            "The {} sampling is missing {} directions."
            .format(state_name, len(missing)))

    # Every direction independently reproduces the same probing base state.
    # Validate its effective stress before assembling the surface.
    base_stresses = np.asarray([
        [float(records[index]["base-sigma-" + name])
         for name in COMPONENT_NAMES]
        for index in range(len(directions))])
    if not np.allclose(
            base_stresses, base_stresses[0], rtol=1.0e-6, atol=1.0e-9):
        raise RuntimeError(
            "The independently reproduced base states do not have "
            "consistent effective stresses.")

    output = os.path.join(root_directory, "yield_surface.csv")
    with open(output, "w", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=CSV_HEADER)
        writer.writeheader()
        for index in range(len(directions)):
            writer.writerow(records[index])

    # Store the sampling grid alongside the ordered surface values so the
    # triangulated surface can be reconstructed without regenerating it.
    np.savetxt(os.path.join(root_directory, "directions.csv"),
               np.column_stack((np.arange(len(directions)), directions)),
               delimiter=",", header="direction-index,n1,n2,n3",
               comments="", fmt=["%d", "%.17e", "%.17e", "%.17e"])
    np.savetxt(os.path.join(root_directory, "faces.csv"), faces,
               delimiter=",", header="vertex-1,vertex-2,vertex-3",
               comments="", fmt="%d")
    np.savetxt(os.path.join(root_directory, "base_stress.csv"),
               base_stresses[:1], delimiter=",")
    print("Wrote {} complete directions and {} faces to {}.".format(
        len(directions), len(faces), output))


def main():
    """Run one sampling worker or merge all completed worker files."""

    # Analysis configuration. Workers are numbered from zero. After every
    # worker has completed, set merge_only to True and run the script once
    # more with the same state name and worker count.
    state_name = "initial"
    worker_index = 0
    number_of_workers = 1
    merge_only = False
    subdivisions = 3
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

    state_name = state_name.lower()
    if state_name not in ("initial", "evolved"):
        raise ValueError("The 3D yield state must be initial or evolved.")
    if number_of_workers < 1:
        raise ValueError("At least one sampling worker is required.")
    if worker_index < 0 or worker_index >= number_of_workers:
        raise ValueError(
            "The worker index must lie in [0, number_of_workers).")
    definition = state_definition(state_name)

    # Sampling requires the calibration generated by the planar-section
    # script. Merge-only runs operate solely on completed worker files.
    if not merge_only:
        calibration_directory = os.path.join(
            here, "yield_surface", "results")
        stiffness_path = os.path.join(
            calibration_directory, "effective_stiffness.csv")
        compliance_path = os.path.join(
            calibration_directory, "effective_compliance.csv")
        if not (os.path.isfile(stiffness_path)
                and os.path.isfile(compliance_path)):
            raise RuntimeError(
                "The effective elastic calibration is unavailable. Run "
                "yield_surface_sections.py first.")
        validate_result_configuration(
            os.path.join(
                calibration_directory, "calibration_configuration.csv"),
            {
                "analysis": "elastic_calibration",
                "elastic_strain": elastic_strain,
                "materials": parameters,
                "mesh": mesh_parameters,
            },
            (stiffness_path, compliance_path))
        compliance = elastic_calibration(
            None, parameters, elastic_strain, calibration_directory)

    directions, faces = icosphere(subdivisions)
    root_directory = os.path.join(
        here, "yield_surface", "surface_3d", state_name)
    part_directory = os.path.join(root_directory, "parts")
    os.makedirs(part_directory, exist_ok=True)
    existing_results = [
        os.path.join(root_directory, "yield_surface.csv"),
        os.path.join(root_directory, "directions.csv"),
        os.path.join(root_directory, "faces.csv"),
        os.path.join(root_directory, "base_strain.csv"),
        os.path.join(root_directory, "base_stress.csv"),
    ]
    existing_results.extend(
        os.path.join(part_directory, file_name)
        for file_name in os.listdir(part_directory)
        if file_name.endswith(".csv"))
    validate_result_configuration(
        os.path.join(root_directory, "sampling_configuration.csv"),
        {
            "analysis": "yield_surface_3d",
            "criterion": criterion,
            "elastic_strain": elastic_strain,
            "materials": parameters,
            "maximum_strain_increment": maximum_strain_increment,
            "mesh": mesh_parameters,
            "number_of_workers": number_of_workers,
            "preload_strain": definition["preload_strain"],
            "base_strain": definition["base_strain"],
            "state": state_name,
            "strain_limits": definition["strain_limits"],
            "subdivisions": subdivisions,
        },
        existing_results)
    if merge_only:
        merge_worker_results(
            state_name, root_directory, directions, faces,
            number_of_workers)
        np.savetxt(os.path.join(root_directory, "base_strain.csv"),
                   definition["base_strain"][None, :], delimiter=",")
        return

    analysis_directory = os.path.join(
        root_directory, "abq_worker_{:02d}".format(worker_index))
    os.makedirs(analysis_directory, exist_ok=True)

    part_path = os.path.join(
        part_directory, "worker_{:02d}.csv".format(worker_index))
    completed = read_worker_results(part_path)
    selected = list(range(
        worker_index, len(directions), number_of_workers))
    invalid = sorted(set(completed).difference(selected))
    if invalid:
        raise RuntimeError(
            "Worker {} contains directions assigned to another worker: {}."
            .format(worker_index, invalid))
    validate_worker_directions(completed, directions)

    original_directory = os.getcwd()
    try:
        os.chdir(analysis_directory)
        gmsh_cell = build_gmsh_cell(
            mesh_parameters,
            "EXAMPLE_3_YIELD3D_{}_{:02d}".format(
                state_name.upper(), worker_index))

        completed_count = len(set(selected).intersection(completed))
        print("{} worker {}/{}: {} assigned directions, {} complete."
              .format(state_name, worker_index + 1, number_of_workers,
                      len(selected), completed_count))
        sys.stdout.flush()

        for position, direction_index in enumerate(selected, start=1):
            if direction_index in completed:
                continue

            # Map the unit normal-stress direction to a proportional strain
            # direction through the effective elastic compliance.
            direction = directions[direction_index]
            stress_direction = np.zeros(6)
            stress_direction[:3] = direction
            strain_direction = compliance.dot(stress_direction)
            strain_direction /= np.max(np.abs(strain_direction))

            print("  worker {}: {}/{} direction {} ({:+.4f}, {:+.4f}, "
                  "{:+.4f})".format(
                      worker_index, position, len(selected), direction_index,
                      direction[0], direction[1], direction[2]))
            sys.stdout.flush()

            for attempt in range(2):
                try:
                    state = sample_direction(
                        gmsh_cell, parameters,
                        definition["preload_strain"],
                        definition["base_strain"], strain_direction,
                        definition["strain_limits"], criterion,
                        maximum_strain_increment)
                    break
                except Exception as error:
                    if attempt == 1:
                        raise
                    print("    analysis failed once; retrying: {}"
                          .format(error))
                    sys.stdout.flush()

            odb_path = state.pop("odb_path")
            append_worker_result(
                part_path, result_row(direction_index, direction, state))
            remove_output_database(odb_path)

        print("*** {} worker {} finished ***".format(
            state_name, worker_index))
    finally:
        os.chdir(original_directory)
        if gmsh.is_initialized():
            gmsh.finalize()


if __name__ == "__main__":
    main()
