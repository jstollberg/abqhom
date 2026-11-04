# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-or-later

"""Verify beam clearance throughout the three Octet load histories.

Every retained beam element in the central cell is compared with the central
cell and its 26 periodic neighbors. Only pairs that overlap within an
idealized Octet joint in the reference configuration are excluded. Contact is
absent while all remaining centerline distances exceed the beam diameter.
This optional post-processing script requires the ``vtk`` Python package and
runs outside Abaqus.
"""

import os
import xml.etree.ElementTree as ET

import numpy as np
import vtk
from octet import octet_topology
from vtk.util.numpy_support import vtk_to_numpy

here = os.getcwd()


def segment_distances(first_start, first_end, second_start, second_end):
    """Return exact distances between corresponding line segments."""
    # Solve the closest-point problem on the supporting lines and clamp both
    # parameters to cover endpoint as well as interior closest points.
    first_direction = first_end - first_start
    second_direction = second_end - second_start
    relative = first_start - second_start

    first_norm = np.einsum("ij,ij->i", first_direction, first_direction)
    second_norm = np.einsum("ij,ij->i", second_direction, second_direction)
    coupling = np.einsum("ij,ij->i", first_direction, second_direction)
    first_projection = np.einsum(
        "ij,ij->i", first_direction, relative)
    second_projection = np.einsum(
        "ij,ij->i", second_direction, relative)
    denominator = first_norm * second_norm - coupling * coupling

    first_parameter = np.zeros_like(denominator)
    regular = denominator > 1.0e-24
    first_parameter[regular] = np.clip(
        (coupling[regular] * second_projection[regular]
         - first_projection[regular] * second_norm[regular])
        / denominator[regular], 0.0, 1.0)
    second_parameter = (
        coupling * first_parameter + second_projection) / second_norm

    below = second_parameter < 0.0
    second_parameter[below] = 0.0
    first_parameter[below] = np.clip(
        -first_projection[below] / first_norm[below], 0.0, 1.0)

    above = second_parameter > 1.0
    second_parameter[above] = 1.0
    first_parameter[above] = np.clip(
        (coupling[above] - first_projection[above]) / first_norm[above],
        0.0, 1.0)

    separation = (relative
                  + first_parameter[:, None] * first_direction
                  - second_parameter[:, None] * second_direction)
    return np.linalg.norm(separation, axis=1)


def read_collection(vtk_root, case_name, model_name):
    """Return ordered pseudo-time/file pairs from a VTK collection."""
    directory = os.path.join(vtk_root, case_name, model_name)
    collection = os.path.join(directory, model_name + ".pvd")
    if not os.path.isfile(collection):
        raise IOError("VTK collection '{}' does not exist.".format(
            collection))

    records = []
    root = ET.parse(collection).getroot()
    for dataset in root.iter("DataSet"):
        records.append((
            float(dataset.attrib["timestep"]),
            os.path.join(directory, dataset.attrib["file"])))
    if not records:
        raise RuntimeError(
            "VTK collection '{}' contains no frames.".format(collection))
    return sorted(records, key=lambda record: record[0])


def read_frame(path):
    """Return reference points, displacements, and line connectivity."""
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(path)
    reader.Update()
    grid = reader.GetOutput()
    if grid.GetPoints() is None:
        raise RuntimeError("VTK frame '{}' contains no points.".format(path))

    displacement_array = grid.GetPointData().GetArray("U")
    if displacement_array is None:
        raise RuntimeError(
            "VTK frame '{}' contains no displacement field 'U'."
            .format(path))

    points = vtk_to_numpy(grid.GetPoints().GetData()).astype(float)
    displacement = vtk_to_numpy(displacement_array).astype(float)
    connectivity = np.empty((grid.GetNumberOfCells(), 2), dtype=int)
    for index in range(grid.GetNumberOfCells()):
        cell = grid.GetCell(index)
        if cell.GetNumberOfPoints() != 2:
            raise RuntimeError(
                "VTK cell {} is not a two-node beam element."
                .format(index))
        connectivity[index] = [cell.GetPointId(0), cell.GetPointId(1)]
    return points, displacement, connectivity


def prescribed_state(case_name, pseudo_time):
    """Return the macroscopic deformation gradient and loading magnitude."""
    paths = {
        "compression": ((0, 0), -0.40),
        "tension": ((0, 0), 0.40),
        "shear_12": ((0, 1), 0.40),
    }
    if case_name not in paths:
        raise ValueError("Unknown load case '{}'.".format(case_name))

    component, final_value = paths[case_name]
    gradient = np.zeros((3, 3))
    gradient[component] = final_value * pseudo_time
    return np.eye(3) + gradient, abs(final_value) * pseudo_time


def assign_elements_to_struts(reference_start, reference_end, edge_length):
    """Assign every imperfect beam element to its perfect parent strut."""
    points, connectivity = octet_topology(edge_length)
    strut_start = points[connectivity[:, 0]]
    strut_end = points[connectivity[:, 1]]
    midpoint = 0.5 * (reference_start + reference_end)
    direction = strut_end - strut_start
    length_squared = np.einsum("ij,ij->i", direction, direction)
    relative = midpoint[:, None, :] - strut_start[None, :, :]
    parameter = (np.einsum("esj,sj->es", relative, direction)
                 / length_squared[None, :])
    parameter = np.clip(parameter, 0.0, 1.0)
    projection = (strut_start[None, :, :]
                  + parameter[:, :, None] * direction[None, :, :])
    distance = np.linalg.norm(midpoint[:, None, :] - projection, axis=2)
    assignment = np.argmin(distance, axis=1)
    maximum_error = np.max(distance[np.arange(len(midpoint)), assignment])
    if maximum_error > 2.0e-3 * edge_length:
        raise RuntimeError(
            "Could not identify all parent struts; maximum error = {}."
            .format(maximum_error))
    return assignment, strut_start, strut_end


def build_admissible_pair_masks(
        reference_start, reference_end, element_strut, strut_start,
        strut_end, integer_shifts, edge_length, diameter):
    """Exclude only reference overlap inherited from ideal joints."""
    element_start = strut_start[element_strut]
    element_end = strut_end[element_strut]
    number_of_elements = len(element_strut)
    masks = np.ones(
        (len(integer_shifts), number_of_elements, number_of_elements),
        dtype=bool)

    reference_midpoint = 0.5 * (reference_start + reference_end)
    reference_length = np.linalg.norm(
        reference_end - reference_start, axis=1)
    half_length_sum = 0.5 * (
        reference_length[:, None] + reference_length[None, :])
    topology_tolerance = 1.0e-8 * edge_length
    contact_tolerance = 1.0e-10 * edge_length

    # Project the imperfect element endpoints onto their perfect parent
    # struts. These projected segments distinguish the overlap inherent to
    # an ideal joint from contact already introduced by the imperfection.
    parent_direction = element_end - element_start
    parent_length_squared = np.einsum(
        "ij,ij->i", parent_direction, parent_direction)
    if np.any(parent_length_squared <= 0.0):
        raise RuntimeError("The Octet topology contains a zero-length strut.")

    start_parameter = np.einsum(
        "ij,ij->i", reference_start - element_start, parent_direction)
    start_parameter /= parent_length_squared
    end_parameter = np.einsum(
        "ij,ij->i", reference_end - element_start, parent_direction)
    end_parameter /= parent_length_squared
    start_parameter = np.clip(start_parameter, 0.0, 1.0)
    end_parameter = np.clip(end_parameter, 0.0, 1.0)
    ideal_start = (
        element_start + start_parameter[:, None] * parent_direction)
    ideal_end = element_start + end_parameter[:, None] * parent_direction
    ideal_length = np.linalg.norm(ideal_end - ideal_start, axis=1)
    if np.any(ideal_length <= topology_tolerance):
        raise RuntimeError(
            "An imperfect beam element does not span a finite interval on "
            "its assigned perfect parent strut.")

    for shift_index, integer_shift in enumerate(integer_shifts):
        image_shift = edge_length * integer_shift
        image_parent_start = element_start + image_shift
        image_parent_end = element_end + image_shift
        junction_distance = np.minimum.reduce((
            np.linalg.norm(
                element_start[:, None, :]
                - image_parent_start[None, :, :],
                axis=2),
            np.linalg.norm(
                element_start[:, None, :]
                - image_parent_end[None, :, :],
                axis=2),
            np.linalg.norm(
                element_end[:, None, :]
                - image_parent_start[None, :, :],
                axis=2),
            np.linalg.norm(
                element_end[:, None, :]
                - image_parent_end[None, :, :],
                axis=2)))

        connected = junction_distance <= topology_tolerance

        # Exact distances are needed only where the midpoint bound permits
        # reference overlap of the circular beam sections.
        image_start = reference_start + image_shift
        image_end = reference_end + image_shift
        image_midpoint = reference_midpoint + image_shift
        midpoint_distance = np.linalg.norm(
            reference_midpoint[:, None, :]
            - image_midpoint[None, :, :], axis=2)
        possible_overlap = (
            midpoint_distance - half_length_sum
            <= diameter + contact_tolerance)
        first, second = np.nonzero(possible_overlap)
        if len(first) == 0:
            continue

        distance = segment_distances(
            reference_start[first], reference_end[first],
            image_start[second], image_end[second])
        overlap = distance <= diameter + contact_tolerance
        if not np.any(overlap):
            continue

        first = first[overlap]
        second = second[overlap]
        overlap_distance = distance[overlap]
        ideal_distance = segment_distances(
            ideal_start[first], ideal_end[first],
            ideal_start[second] + image_shift,
            ideal_end[second] + image_shift)
        intentional = (connected[first, second]
                       & (ideal_distance
                          <= diameter + contact_tolerance))
        if np.any(np.invert(intentional)):
            local = int(np.flatnonzero(np.invert(intentional))[0])
            raise RuntimeError(
                "Beam elements overlap outside an intentional ideal-joint "
                "region in the reference cell or one of its periodic "
                "images: elements ({}, {}), periodic shift {}, centerline "
                "distance {:.8e}, diameter {:.8e}.".format(
                    int(first[local]), int(second[local]), integer_shift,
                    overlap_distance[local], diameter))

        # Exempt only the element pairs that already form the finite overlap
        # region of an ideal rigid joint. Portions outside this reference
        # overlap remain contact candidates even on connected struts.
        masks[shift_index, first, second] = False

    return masks


def minimum_periodic_distance(current_start, current_end, current_shifts,
                              admissible_pairs, seed_pair=None):
    """Return the exact minimum distance over all admissible image pairs."""
    midpoint = 0.5 * (current_start + current_end)
    length = np.linalg.norm(current_end - current_start, axis=1)
    half_length_sum = 0.5 * (length[:, None] + length[None, :])

    minimum = np.inf
    critical_pair = None

    # The previous frame's critical pair supplies an inexpensive finite
    # upper bound. It is not used as a lower bound and therefore cannot
    # exclude a newly closer pair.
    if seed_pair is not None:
        first, shift_index, second = seed_pair
        if admissible_pairs[shift_index, first, second]:
            image_shift = current_shifts[shift_index]
            minimum = segment_distances(
                current_start[first][None, :],
                current_end[first][None, :],
                (current_start[second] + image_shift)[None, :],
                (current_end[second] + image_shift)[None, :])[0]
            critical_pair = seed_pair

    for shift_index, image_shift in enumerate(current_shifts):
        image_start = current_start + image_shift
        image_end = current_end + image_shift
        image_midpoint = midpoint + image_shift

        midpoint_distance = np.linalg.norm(
            midpoint[:, None, :] - image_midpoint[None, :, :], axis=2)

        # Every point of a segment lies within half its length of its
        # midpoint. Hence this is a rigorous lower bound on segment distance.
        lower_bound = midpoint_distance - half_length_sum

        # Only pairs whose lower bound can improve the current upper bound
        # need the more expensive exact segment-distance evaluation.
        candidates = (admissible_pairs[shift_index]
                      & (lower_bound < minimum))
        first, second = np.nonzero(candidates)
        if len(first) == 0:
            continue

        distance = segment_distances(
            current_start[first], current_end[first],
            image_start[second], image_end[second])
        local = np.argmin(distance)
        if distance[local] < minimum:
            minimum = distance[local]
            critical_pair = (
                int(first[local]), shift_index, int(second[local]))

    if critical_pair is None:
        raise RuntimeError("No nonconnected beam-element pair was found.")
    return minimum, critical_pair


def evaluate_clearance(
        case_name, vtk_root, model_name, edge_length, diameter):
    """Evaluate all non-overlapping centerline pairs throughout one case."""
    records = read_collection(vtk_root, case_name, model_name)
    reference_points, _, connectivity = read_frame(records[0][1])
    reference_start = reference_points[connectivity[:, 0]]
    reference_end = reference_points[connectivity[:, 1]]
    element_strut, strut_start, strut_end = assign_elements_to_struts(
        reference_start, reference_end, edge_length)

    integer_shifts = np.asarray(
        [[0, 0, 0]] + [
            [i, j, k]
            for i in (-1, 0, 1)
            for j in (-1, 0, 1)
            for k in (-1, 0, 1)
            if (i, j, k) != (0, 0, 0)], dtype=float)
    admissible_pairs = build_admissible_pair_masks(
        reference_start, reference_end, element_strut, strut_start,
        strut_end, integer_shifts, edge_length, diameter)

    output = []
    global_minimum = np.inf
    closest_pair = None
    previous_pair = None
    for frame_index, (pseudo_time, path) in enumerate(records):
        points, displacement, frame_connectivity = read_frame(path)
        if not np.array_equal(frame_connectivity, connectivity):
            raise RuntimeError("VTK connectivity changed between frames.")

        current_points = points + displacement
        current_start = current_points[connectivity[:, 0]]
        current_end = current_points[connectivity[:, 1]]

        deformation, loading_magnitude = prescribed_state(
            case_name, pseudo_time)

        # Periodic fluctuations repeat, while each reference lattice vector
        # is carried to the current configuration by the macroscopic F.
        current_shifts = edge_length * integer_shifts.dot(deformation.T)
        minimum, frame_pair = minimum_periodic_distance(
            current_start, current_end, current_shifts,
            admissible_pairs, seed_pair=previous_pair)
        previous_pair = frame_pair

        if minimum < global_minimum:
            global_minimum = minimum
            first, shift_index, second = frame_pair
            image_shift = edge_length * integer_shifts[shift_index]
            reference_distance = segment_distances(
                reference_start[first][None, :],
                reference_end[first][None, :],
                (reference_start[second] + image_shift)[None, :],
                (reference_end[second] + image_shift)[None, :])[0]
            closest_pair = {
                "frame": frame_index,
                "pseudo_time": pseudo_time,
                "first_cell": first,
                "second_cell": second,
                "shift": integer_shifts[shift_index].copy(),
                "reference_distance": reference_distance,
            }
        output.append([pseudo_time, loading_magnitude, minimum])

        if ((frame_index + 1) % 50 == 0
                or frame_index + 1 == len(records)):
            print("{}: checked frame {} of {}".format(
                case_name, frame_index + 1, len(records)))

    return np.asarray(output), closest_pair


def report_case(case_name, result, closest_pair, diameter):
    """Print the limiting clearance and corresponding element pair."""
    frame = np.argmin(result[:, 2])
    distance = result[frame, 2]
    print("\n{}:".format(case_name))
    print("  minimum centerline distance: {:.8e}".format(distance))
    print("  minimum surface clearance  : {:.8e}".format(
        distance - diameter))
    print("  loading magnitude          : {:.6f}".format(
        result[frame, 1]))

    if closest_pair is None:
        print("  no independent pair entered the near-contact search radius")
        return

    print("  critical frame             : {} (tau = {:.6f})".format(
        closest_pair["frame"], closest_pair["pseudo_time"]))
    print("  VTK cell indices (zero-based): ({}, {})".format(
        closest_pair["first_cell"], closest_pair["second_cell"]))
    print("  periodic image shift       : {}".format(
        closest_pair["shift"]))
    print("  initial centerline distance: {:.8e}".format(
        closest_pair["reference_distance"]))


def assert_no_contact(case_name, result, closest_pair, diameter,
                      edge_length):
    """Raise an error if two nonconnected beam surfaces touch or overlap."""
    frame = np.argmin(result[:, 2])
    distance = result[frame, 2]
    tolerance = 1.0e-10 * edge_length
    if distance > diameter + tolerance:
        return

    # Tangency and positive overlap both invalidate the no-contact model;
    # values within tolerance are treated conservatively as unresolved.
    raise RuntimeError(
        "Contact detected or unresolved in '{}' at tau={:.8e}: centerline "
        "distance {:.8e}, diameter {:.8e}, elements ({}, {}), periodic "
        "shift {}.".format(
            case_name, closest_pair["pseudo_time"], distance, diameter,
            closest_pair["first_cell"], closest_pair["second_cell"],
            closest_pair["shift"]))


def main():
    """Check all stored load histories and export clearance curves."""
    edge_length = 1.0
    beam_radius = 0.04
    diameter = 2.0 * beam_radius
    model_name = "EXAMPLE_2_OCTET"
    cases = ("compression", "tension", "shear_12")
    vtk_root = os.path.join(here, "vtk_output")
    result_root = os.path.join(here, "results")

    for case_name in cases:
        result, closest_pair = evaluate_clearance(
            case_name, vtk_root, model_name, edge_length, diameter)
        result_directory = os.path.join(result_root, case_name)
        if not os.path.isdir(result_directory):
            os.makedirs(result_directory)
        np.savetxt(
            os.path.join(result_directory, "contact_clearance.csv"),
            result, delimiter=",")
        report_case(case_name, result, closest_pair, diameter)
        assert_no_contact(
            case_name, result, closest_pair, diameter, edge_length)


if __name__ == "__main__":
    main()
