# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-or-later

"""Verify beam clearance throughout the three Octet load histories.

Every retained beam element in the central cell is compared with the central
cell and its 26 periodic neighbors. Centerline intervals of connected members
inside a geometry-based spherical neighborhood of their reference idealized
joint are excluded. These material intervals remain excluded throughout the
history. All other intervals remain contact candidates. Contact is absent
while all remaining centerline distances exceed the beam diameter.
This optional post-processing script requires the ``vtk`` Python package and
runs outside Abaqus.
"""

import os
import xml.etree.ElementTree as ET

import numpy as np
import vtk
from octet import octet_topology
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

from abqhom.utils import export_to_csv_file

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


def segment_intervals_outside_sphere(
        start, end, center, radius, tolerance):
    """Return segment-parameter intervals outside a reference sphere."""
    direction = end - start
    length_squared = np.dot(direction, direction)
    if length_squared <= tolerance * tolerance:
        raise RuntimeError("A zero-length beam segment was encountered.")

    relative = start - center
    linear = 2.0 * np.dot(relative, direction)
    constant = np.dot(relative, relative) - radius * radius
    discriminant = linear * linear - 4.0 * length_squared * constant

    parameters = [0.0, 1.0]
    if discriminant > tolerance * tolerance:
        root = np.sqrt(discriminant)
        for value in ((-linear - root) / (2.0 * length_squared),
                      (-linear + root) / (2.0 * length_squared)):
            if 0.0 < value < 1.0:
                parameters.append(value)
    parameters = sorted(parameters)

    intervals = []
    parameter_tolerance = tolerance / np.sqrt(length_squared)
    for lower, upper in zip(parameters[:-1], parameters[1:]):
        if upper - lower <= parameter_tolerance:
            continue
        midpoint = start + 0.5 * (lower + upper) * direction
        if np.linalg.norm(midpoint - center) <= radius + tolerance:
            continue
        intervals.append((lower, upper))
    return intervals


def connected_segment_distance(first_start, first_end, second_start,
                               second_end, first_intervals,
                               second_intervals):
    """Return segment distance outside their common idealized joint."""
    if not first_intervals or not second_intervals:
        return np.inf

    first_direction = first_end - first_start
    second_direction = second_end - second_start
    first = []
    second = []
    for first_lower, first_upper in first_intervals:
        first_portion = (
            first_start + first_lower * first_direction,
            first_start + first_upper * first_direction)
        for second_lower, second_upper in second_intervals:
            second_portion = (
                second_start + second_lower * second_direction,
                second_start + second_upper * second_direction)
            first.append(first_portion)
            second.append(second_portion)
    first = np.asarray(first)
    second = np.asarray(second)
    return np.min(segment_distances(
        first[:, 0], first[:, 1], second[:, 0], second[:, 1]))


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


def octet_joint_radius(edge_length, diameter):
    """Return a joint radius enclosing nominal connected-beam overlap."""
    points, connectivity = octet_topology(edge_length)
    maximum_extent = diameter
    for point_index, point in enumerate(points):
        directions = []
        for first, second in connectivity:
            if first == point_index:
                directions.append(points[second] - point)
            elif second == point_index:
                directions.append(points[first] - point)
        directions = np.asarray(directions, dtype=float)
        directions /= np.linalg.norm(directions, axis=1)[:, None]

        for first in range(len(directions)):
            for second in range(first + 1, len(directions)):
                cosine = np.clip(
                    np.dot(directions[first], directions[second]),
                    -1.0, 1.0)
                if cosine >= 1.0 - 1.0e-12:
                    raise RuntimeError(
                        "The Octet topology contains coincident incident "
                        "struts at one joint.")
                if cosine <= -1.0 + 1.0e-12:
                    # Oppositely directed collinear struts form a continuous
                    # member rather than a finite-angle joint overlap.
                    continue
                angle = np.arccos(cosine)
                effective_angle = min(angle, np.pi - angle)
                sine = np.sin(effective_angle)
                maximum_extent = max(maximum_extent, diameter / sine)
    return maximum_extent


def map_strut_endpoints(reference_start, reference_end, element_strut,
                        strut_start, strut_end, edge_length):
    """Map active parent-strut endpoints to element endpoint indices."""
    number_of_struts = len(strut_start)
    element_index = -np.ones((number_of_struts, 2), dtype=int)
    endpoint_index = -np.ones((number_of_struts, 2), dtype=int)
    element_points = np.stack((reference_start, reference_end), axis=1)
    parent_points = np.stack((strut_start, strut_end), axis=1)
    tolerance = 2.0e-3 * edge_length

    for strut in np.unique(element_strut):
        candidates = np.flatnonzero(element_strut == strut)
        candidate_points = element_points[candidates]
        for endpoint in (0, 1):
            distance = np.linalg.norm(
                candidate_points - parent_points[strut, endpoint], axis=2)
            local = int(np.argmin(distance))
            candidate, side = np.unravel_index(local, distance.shape)
            if distance[candidate, side] > tolerance:
                raise RuntimeError(
                    "Could not locate an active endpoint of parent strut "
                    "{}.".format(strut))
            element_index[strut, endpoint] = candidates[candidate]
            endpoint_index[strut, endpoint] = side
    return element_index, endpoint_index


def strut_endpoint_coordinates(start, end, element_index, endpoint_index):
    """Return endpoint coordinates of all active parent struts."""
    endpoints = np.full(element_index.shape + (3,), np.nan)
    for strut in range(len(element_index)):
        for endpoint in (0, 1):
            element = element_index[strut, endpoint]
            side = endpoint_index[strut, endpoint]
            if element < 0:
                continue
            endpoints[strut, endpoint] = (
                start[element] if side == 0 else end[element])
    return endpoints


def build_joint_exterior_intervals(
        reference_start, reference_end, element_strut,
        reference_strut_endpoints, joint_radius, tolerance):
    """Identify each element's material portions outside either joint."""
    intervals = []
    for element, strut in enumerate(element_strut):
        element_intervals = []
        for endpoint in (0, 1):
            center = reference_strut_endpoints[strut, endpoint]
            if not np.all(np.isfinite(center)):
                raise RuntimeError(
                    "Could not locate a reference idealized joint.")
            element_intervals.append(segment_intervals_outside_sphere(
                reference_start[element], reference_end[element], center,
                joint_radius, tolerance))
        intervals.append(element_intervals)
    return intervals


def build_pair_topology(
        reference_start, reference_end, element_strut, strut_start,
        strut_end, integer_shifts, edge_length, diameter):
    """Build masks and common-joint endpoints for periodic segment pairs."""
    number_of_elements = len(element_strut)
    number_of_struts = len(strut_start)
    masks = np.ones(
        (len(integer_shifts), number_of_elements, number_of_elements),
        dtype=bool)
    first_joint = -np.ones(masks.shape, dtype=np.int8)
    second_joint = -np.ones(masks.shape, dtype=np.int8)
    parent_points = np.stack((strut_start, strut_end), axis=1)
    topology_tolerance = 1.0e-8 * edge_length
    contact_tolerance = 1.0e-10 * edge_length

    first_strut = element_strut[:, None]
    second_strut = element_strut[None, :]
    for shift_index, integer_shift in enumerate(integer_shifts):
        image_shift = edge_length * integer_shift
        image_points = parent_points + image_shift
        endpoint_distance = np.linalg.norm(
            parent_points[:, None, :, None, :]
            - image_points[None, :, None, :, :], axis=4)
        matching = endpoint_distance <= topology_tolerance
        number_of_matches = np.sum(matching, axis=(2, 3))

        # One matching endpoint identifies a common idealized joint. Two
        # matches identify the same parent strut or a periodic image of it.
        connected = number_of_matches == 1
        flattened = np.argmax(matching.reshape(
            number_of_struts, number_of_struts, 4), axis=2)
        parent_first_joint = np.where(connected, flattened // 2, -1)
        parent_second_joint = np.where(connected, flattened % 2, -1)
        first_joint[shift_index] = parent_first_joint[
            first_strut, second_strut]
        second_joint[shift_index] = parent_second_joint[
            first_strut, second_strut]

        coincident = (number_of_matches >= 2)[first_strut, second_strut]
        first, second = np.nonzero(coincident)
        if len(first) == 0:
            continue

        distance = segment_distances(
            reference_start[first], reference_end[first],
            reference_start[second] + image_shift,
            reference_end[second] + image_shift)
        local_overlap = distance <= diameter + contact_tolerance
        masks[shift_index, first[local_overlap],
              second[local_overlap]] = False

    return masks, first_joint, second_joint


def pair_distance(first, shift_index, second, current_start, current_end,
                  current_shifts, first_joint, second_joint,
                  exterior_intervals):
    """Return one periodic segment-pair distance outside ideal joints."""
    image_shift = current_shifts[shift_index]
    second_start = current_start[second] + image_shift
    second_end = current_end[second] + image_shift
    first_endpoint = first_joint[shift_index, first, second]
    if first_endpoint < 0:
        return segment_distances(
            current_start[first][None, :],
            current_end[first][None, :],
            second_start[None, :], second_end[None, :])[0]

    second_endpoint = second_joint[shift_index, first, second]
    return connected_segment_distance(
        current_start[first], current_end[first], second_start, second_end,
        exterior_intervals[first][first_endpoint],
        exterior_intervals[second][second_endpoint])


def minimum_periodic_distance(
        current_start, current_end, current_shifts, admissible_pairs,
        first_joint, second_joint, exterior_intervals, seed_pair=None):
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
            minimum = pair_distance(
                first, shift_index, second, current_start, current_end,
                current_shifts, first_joint, second_joint,
                exterior_intervals)
            if np.isfinite(minimum):
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

        connected = first_joint[shift_index, first, second] >= 0

        # Nonconnected segments need no clipping and remain vectorized.
        independent = np.invert(connected)
        if np.any(independent):
            independent_first = first[independent]
            independent_second = second[independent]
            distance = segment_distances(
                current_start[independent_first],
                current_end[independent_first],
                image_start[independent_second],
                image_end[independent_second])
            local = np.argmin(distance)
            if distance[local] < minimum:
                minimum = distance[local]
                critical_pair = (
                    int(independent_first[local]), shift_index,
                    int(independent_second[local]))

        if not np.any(connected):
            continue

        # The unclipped distance is a lower bound for the distance after
        # removing the joint portions. It cheaply rejects connected pairs
        # that cannot improve the current minimum.
        connected_first = first[connected]
        connected_second = second[connected]
        lower = segment_distances(
            current_start[connected_first], current_end[connected_first],
            image_start[connected_second], image_end[connected_second])
        possible = np.flatnonzero(lower < minimum)
        for local in possible:
            first_element = int(connected_first[local])
            second_element = int(connected_second[local])
            distance = pair_distance(
                first_element, shift_index, second_element,
                current_start, current_end, current_shifts, first_joint,
                second_joint, exterior_intervals)
            if distance < minimum:
                minimum = distance
                critical_pair = (
                    first_element, shift_index, second_element)

    if critical_pair is None:
        raise RuntimeError(
            "No beam-segment pair outside the idealized joints was found.")
    return minimum, critical_pair


def find_contact_pairs(
        current_start, current_end, current_shifts, admissible_pairs,
        first_joint, second_joint, exterior_intervals, diameter,
        tolerance):
    """Return all exterior segment pairs within the contact threshold."""
    midpoint = 0.5 * (current_start + current_end)
    length = np.linalg.norm(current_end - current_start, axis=1)
    half_length_sum = 0.5 * (length[:, None] + length[None, :])
    threshold = diameter + tolerance
    contacts = []

    for shift_index, image_shift in enumerate(current_shifts):
        image_start = current_start + image_shift
        image_end = current_end + image_shift
        image_midpoint = midpoint + image_shift
        midpoint_distance = np.linalg.norm(
            midpoint[:, None, :] - image_midpoint[None, :, :], axis=2)
        lower_bound = midpoint_distance - half_length_sum
        first, second = np.nonzero(
            admissible_pairs[shift_index] & (lower_bound <= threshold))
        if len(first) == 0:
            continue

        connected = first_joint[shift_index, first, second] >= 0
        independent = np.invert(connected)
        if np.any(independent):
            independent_first = first[independent]
            independent_second = second[independent]
            distance = segment_distances(
                current_start[independent_first],
                current_end[independent_first],
                image_start[independent_second],
                image_end[independent_second])
            for local in np.flatnonzero(distance <= threshold):
                contacts.append((
                    int(independent_first[local]), shift_index,
                    int(independent_second[local])))

        for local in np.flatnonzero(connected):
            first_element = int(first[local])
            second_element = int(second[local])
            distance = pair_distance(
                first_element, shift_index, second_element,
                current_start, current_end, current_shifts,
                first_joint, second_joint, exterior_intervals)
            if distance <= threshold:
                contacts.append(
                    (first_element, shift_index, second_element))

    return contacts


def export_contact_state(
        path, current_points, connectivity, integer_shifts,
        current_shifts, contacts):
    """Export the critical configuration and mark contacting elements."""
    if not os.path.isdir(os.path.dirname(path)):
        os.makedirs(os.path.dirname(path))

    central_flags = np.zeros(len(connectivity), dtype=np.int32)
    image_flags = {}
    for first, shift_index, second in contacts:
        central_flags[first] = 1
        if shift_index == 0:
            central_flags[second] = 1
        else:
            flags = image_flags.setdefault(
                shift_index, np.zeros(len(connectivity), dtype=np.int32))
            flags[second] = 1

    displayed_shifts = [0] + sorted(image_flags)
    point_blocks = []
    contact_blocks = []
    source_blocks = []
    shift_blocks = []
    cell_blocks = []
    point_offset = 0
    for shift_index in displayed_shifts:
        point_blocks.append(current_points + current_shifts[shift_index])
        cell_blocks.append(connectivity + point_offset)
        point_offset += len(current_points)

        flags = (central_flags if shift_index == 0
                 else image_flags[shift_index])
        contact_blocks.append(flags)
        source_blocks.append(
            np.arange(len(connectivity), dtype=np.int32))
        shift = np.asarray(integer_shifts[shift_index], dtype=np.int32)
        shift_blocks.append(np.tile(shift, (len(connectivity), 1)))

    grid = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()
    points.SetData(numpy_to_vtk(
        np.concatenate(point_blocks), deep=True))
    grid.SetPoints(points)

    for first, second in np.concatenate(cell_blocks):
        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, int(first))
        line.GetPointIds().SetId(1, int(second))
        grid.InsertNextCell(line.GetCellType(), line.GetPointIds())

    def add_cell_data(name, values, active=False):
        """Attach one NumPy array as VTK cell data."""
        array = numpy_to_vtk(values, deep=True)
        array.SetName(name)
        if active:
            grid.GetCellData().SetScalars(array)
        else:
            grid.GetCellData().AddArray(array)

    add_cell_data("CONTACT", np.concatenate(contact_blocks), active=True)
    add_cell_data("SOURCE_ELEMENT", np.concatenate(source_blocks))
    shifts = np.concatenate(shift_blocks)
    add_cell_data("IMAGE_I", shifts[:, 0])
    add_cell_data("IMAGE_J", shifts[:, 1])
    add_cell_data("IMAGE_K", shifts[:, 2])

    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(path)
    writer.SetInputData(grid)
    if writer.Write() != 1:
        raise IOError("Could not write contact VTK file '{}'.".format(path))


def evaluate_clearance(
        case_name, vtk_root, model_name, edge_length, diameter,
        joint_radius, diagnostic_directory):
    """Evaluate exterior centerline distances throughout one load case."""
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
    admissible_pairs, first_joint, second_joint = build_pair_topology(
        reference_start, reference_end, element_strut, strut_start,
        strut_end, integer_shifts, edge_length, diameter)
    endpoint_element, endpoint_side = map_strut_endpoints(
        reference_start, reference_end, element_strut, strut_start,
        strut_end, edge_length)
    reference_strut_endpoints = strut_endpoint_coordinates(
        reference_start, reference_end, endpoint_element, endpoint_side)
    tolerance = 1.0e-10 * edge_length
    exterior_intervals = build_joint_exterior_intervals(
        reference_start, reference_end, element_strut,
        reference_strut_endpoints, joint_radius, tolerance)
    reference_shifts = edge_length * integer_shifts

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
            admissible_pairs, first_joint, second_joint, exterior_intervals,
            seed_pair=previous_pair)
        previous_pair = frame_pair

        if minimum < global_minimum:
            global_minimum = minimum
            first, shift_index, second = frame_pair
            reference_distance = pair_distance(
                first, shift_index, second,
                reference_start, reference_end, reference_shifts,
                first_joint, second_joint, exterior_intervals)
            closest_pair = {
                "frame": frame_index,
                "pseudo_time": pseudo_time,
                "first_element": first,
                "second_element": second,
                "shift": integer_shifts[shift_index].copy(),
                "reference_distance": reference_distance,
            }
        output.append([pseudo_time, loading_magnitude, minimum])

        if ((frame_index + 1) % 50 == 0
                or frame_index + 1 == len(records)):
            print("{}: checked frame {} of {}".format(
                case_name, frame_index + 1, len(records)))

    critical_frame = closest_pair["frame"]
    points, displacement, frame_connectivity = read_frame(
        records[critical_frame][1])
    if not np.array_equal(frame_connectivity, connectivity):
        raise RuntimeError("VTK connectivity changed between frames.")
    current_points = points + displacement
    current_start = current_points[connectivity[:, 0]]
    current_end = current_points[connectivity[:, 1]]
    deformation, _ = prescribed_state(
        case_name, records[critical_frame][0])
    current_shifts = edge_length * integer_shifts.dot(deformation.T)
    contacts = find_contact_pairs(
        current_start, current_end, current_shifts, admissible_pairs,
        first_joint, second_joint, exterior_intervals, diameter,
        tolerance)
    print("{}: export minimum-clearance configuration.".format(case_name))
    export_contact_state(
        os.path.join(diagnostic_directory, case_name + ".vtu"),
        current_points, connectivity, integer_shifts, current_shifts,
        contacts)

    return np.asarray(output), closest_pair


def report_case(case_name, result, closest_pair, diameter):
    """Print the limiting clearance and corresponding element pair."""
    frame = np.argmin(result[:, 2])
    distance = result[frame, 2]
    print("\n{}:".format(case_name))
    print("  minimum exterior centerline distance: {:.8e}".format(distance))
    print("  minimum exterior surface clearance  : {:.8e}".format(
        distance - diameter))
    print("  loading magnitude                  : {:.6f}".format(
        result[frame, 1]))

    if closest_pair is None:
        print("  no exterior pair entered the near-contact search radius")
        return

    print("  critical frame                     : {} (tau = {:.6f})".format(
        closest_pair["frame"], closest_pair["pseudo_time"]))
    print("  VTK cell indices (zero-based)      : ({}, {})".format(
        closest_pair["first_element"], closest_pair["second_element"]))
    print("  periodic image shift               : {}".format(
        closest_pair["shift"]))
    if np.isfinite(closest_pair["reference_distance"]):
        print("  initial exterior centerline distance: {:.8e}".format(
            closest_pair["reference_distance"]))
    else:
        print("  initial exterior centerline distance: not applicable")


def assert_no_contact(case_name, result, closest_pair, diameter,
                      edge_length):
    """Raise if beam surfaces touch outside idealized joint regions."""
    frame = np.argmin(result[:, 2])
    distance = result[frame, 2]
    tolerance = 1.0e-10 * edge_length
    if distance > diameter + tolerance:
        return

    # Tangency and positive overlap both invalidate the no-contact model;
    # values within tolerance are treated conservatively as unresolved.
    raise RuntimeError(
        "Contact detected or unresolved outside an idealized joint in '{}' "
        "at tau={:.8e}: centerline distance {:.8e}, diameter {:.8e}, "
        "elements ({}, {}), periodic shift {}.".format(
            case_name, closest_pair["pseudo_time"], distance, diameter,
            closest_pair["first_element"], closest_pair["second_element"],
            closest_pair["shift"]))


def main():
    """Check all stored load histories and export clearance curves."""
    edge_length = 1.0
    beam_radius = 0.04
    diameter = 2.0 * beam_radius
    joint_radius = octet_joint_radius(edge_length, diameter)
    model_name = "RVE"
    cases = ("compression", "tension", "shear_12")
    vtk_root = os.path.join(here, "vtk_output")
    diagnostic_root = os.path.join(here, "contact_vtk_output")
    result_root = os.path.join(here, "results")
    print("Idealized joint exclusion radius: {:.8e}".format(joint_radius))

    failures = []
    for case_name in cases:
        result, closest_pair = evaluate_clearance(
            case_name, vtk_root, model_name, edge_length, diameter,
            joint_radius, diagnostic_root)
        result_directory = os.path.join(result_root, case_name)
        if not os.path.isdir(result_directory):
            os.makedirs(result_directory)
        export_to_csv_file(
            result,
            os.path.join(result_directory, "contact_clearance.csv"),
            header=("pseudo_time,loading_measure,"
                    "minimum_exterior_centerline_distance"))
        report_case(case_name, result, closest_pair, diameter)
        try:
            assert_no_contact(
                case_name, result, closest_pair, diameter, edge_length)
        except RuntimeError as error:
            failures.append(str(error))

    if failures:
        raise RuntimeError(
            "Contact verification failed:\n  " + "\n  ".join(failures))


if __name__ == "__main__":
    main()
