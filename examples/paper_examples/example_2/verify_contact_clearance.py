"""Check beam-surface clearance over all final Octet load histories.

The check compares every beam element in the central reduced cell with beam
elements in the central cell and its 26 periodic neighbors. Element pairs
whose centerlines were already close in the reference configuration are
excluded as parts of the same idealized lattice junction. For all remaining
pairs, contact would require a centerline distance below twice the beam
radius.
"""

import os
import xml.etree.ElementTree as ET

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


HERE = os.path.dirname(os.path.abspath(__file__))
EDGE_LENGTH = 1.0
BEAM_RADIUS = 0.04
DIAMETER = 2.0*BEAM_RADIUS
STEP_DURATION = 2.0


def segment_distances(first_start, first_end, second_start, second_end):
    """Return exact distances between corresponding line segments."""
    first_direction = first_end - first_start
    second_direction = second_end - second_start
    relative = first_start - second_start

    first_norm = np.einsum("ij,ij->i", first_direction, first_direction)
    second_norm = np.einsum("ij,ij->i", second_direction, second_direction)
    coupling = np.einsum(
        "ij,ij->i", first_direction, second_direction)
    first_projection = np.einsum(
        "ij,ij->i", first_direction, relative)
    second_projection = np.einsum(
        "ij,ij->i", second_direction, relative)
    denominator = first_norm*second_norm - coupling*coupling

    first_parameter = np.zeros_like(denominator)
    regular = denominator > 1.0e-24
    first_parameter[regular] = np.clip(
        ((coupling[regular]*second_projection[regular]
          - first_projection[regular]*second_norm[regular])
         / denominator[regular]), 0.0, 1.0)
    second_parameter = (
        coupling*first_parameter + second_projection)/second_norm

    below = second_parameter < 0.0
    second_parameter[below] = 0.0
    first_parameter[below] = np.clip(
        -first_projection[below]/first_norm[below], 0.0, 1.0)

    above = second_parameter > 1.0
    second_parameter[above] = 1.0
    first_parameter[above] = np.clip(
        ((coupling[above] - first_projection[above])
         / first_norm[above]), 0.0, 1.0)

    separation = (relative
                  + first_parameter[:, None]*first_direction
                  - second_parameter[:, None]*second_direction)
    return np.linalg.norm(separation, axis=1)


def read_collection(case_name):
    """Return ordered time/file pairs from one Paraqus collection."""
    directory = os.path.join(
        HERE, "final_vtk", case_name, "EXAMPLE_2_OCTET")
    collection = os.path.join(directory, "EXAMPLE_2_OCTET.pvd")
    root = ET.parse(collection).getroot()
    records = []
    for dataset in root.iter("DataSet"):
        records.append((
            float(dataset.attrib["timestep"]),
            os.path.join(directory, dataset.attrib["file"])))
    return records


def read_frame(path):
    """Return reference points, displacements, and line connectivity."""
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(path)
    reader.Update()
    grid = reader.GetOutput()
    points = vtk_to_numpy(grid.GetPoints().GetData()).astype(float)
    displacement = vtk_to_numpy(
        grid.GetPointData().GetArray("U")).astype(float)
    connectivity = np.empty((grid.GetNumberOfCells(), 2), dtype=int)
    for index in range(grid.GetNumberOfCells()):
        cell = grid.GetCell(index)
        connectivity[index] = [cell.GetPointId(0), cell.GetPointId(1)]
    return points, displacement, connectivity


def deformation_gradient(case_name, time):
    """Return the prescribed macroscopic deformation gradient."""
    gradient = np.zeros((3, 3))
    magnitude = 0.4*time/STEP_DURATION
    if case_name == "compression":
        gradient[0, 0] = -magnitude
    elif case_name == "tension":
        gradient[0, 0] = magnitude
    elif case_name == "shear_12":
        gradient[0, 1] = magnitude
    else:
        raise ValueError("Unknown case '{}'".format(case_name))
    return np.eye(3) + gradient


def ideal_octet_struts():
    """Return endpoints of the 36 straight struts in the perfect cell."""
    half_edge = 0.5*EDGE_LENGTH
    corners = {}
    for sign_x in (-1, 1):
        for sign_y in (-1, 1):
            for sign_z in (-1, 1):
                signs = (sign_x, sign_y, sign_z)
                corners[signs] = np.asarray(
                    signs, dtype=float)*half_edge
    face_centers = {}
    for axis in range(3):
        for sign in (-1, 1):
            coordinates = np.zeros(3)
            coordinates[axis] = sign*half_edge
            face_centers[(axis, sign)] = coordinates

    struts = []
    for (axis, sign), face_center in face_centers.items():
        for signs, corner in corners.items():
            if signs[axis] == sign:
                struts.append((face_center, corner))
    face_items = sorted(face_centers.items())
    for first_index, ((first_axis, _), first) in enumerate(face_items):
        for (second_axis, _), second in face_items[first_index + 1:]:
            if first_axis != second_axis:
                struts.append((first, second))
    return (np.asarray([strut[0] for strut in struts]),
            np.asarray([strut[1] for strut in struts]))


def assign_elements_to_struts(reference_start, reference_end):
    """Assign each imperfect beam element to its perfect parent strut."""
    strut_start, strut_end = ideal_octet_struts()
    midpoint = 0.5*(reference_start + reference_end)
    direction = strut_end - strut_start
    length_squared = np.einsum("ij,ij->i", direction, direction)
    relative = midpoint[:, None, :] - strut_start[None, :, :]
    parameter = np.einsum(
        "esj,sj->es", relative, direction)/length_squared[None, :]
    parameter = np.clip(parameter, 0.0, 1.0)
    projection = (strut_start[None, :, :]
                  + parameter[:, :, None]*direction[None, :, :])
    distance = np.linalg.norm(midpoint[:, None, :] - projection, axis=2)
    assignment = np.argmin(distance, axis=1)
    maximum_error = np.max(distance[np.arange(len(midpoint)), assignment])
    if maximum_error > 2.0e-3:
        raise RuntimeError(
            "Could not identify all parent struts; maximum error = {}."
            .format(maximum_error))
    return assignment, strut_start, strut_end


def minimum_clearances(case_name):
    """Evaluate minimum non-junction centerline distances for every frame."""
    records = read_collection(case_name)
    reference_points, _, connectivity = read_frame(records[0][1])
    reference_start = reference_points[connectivity[:, 0]]
    reference_end = reference_points[connectivity[:, 1]]
    number_of_elements = len(connectivity)
    element_strut, strut_start, strut_end = assign_elements_to_struts(
        reference_start, reference_end)

    integer_shifts = np.asarray([
        [i, j, k]
        for i in (-1, 0, 1)
        for j in (-1, 0, 1)
        for k in (-1, 0, 1)], dtype=float)
    zero_shift = np.where(np.all(integer_shifts == 0.0, axis=1))[0][0]

    output = []
    global_minimum = np.inf
    closest_pair = None
    for frame_index, (time, path) in enumerate(records):
        points, displacement, frame_connectivity = read_frame(path)
        if not np.array_equal(frame_connectivity, connectivity):
            raise RuntimeError("VTK connectivity changed between frames.")
        current_points = points + displacement
        current_start = current_points[connectivity[:, 0]]
        current_end = current_points[connectivity[:, 1]]
        current_midpoint = 0.5*(current_start + current_end)
        maximum_length = np.max(np.linalg.norm(
            current_end - current_start, axis=1))

        deformation = deformation_gradient(case_name, time)
        current_shifts = EDGE_LENGTH*integer_shifts.dot(deformation.T)
        copied_start = (
            current_start[None, :, :] + current_shifts[:, None, :])
        copied_end = (
            current_end[None, :, :] + current_shifts[:, None, :])
        copied_midpoint = 0.5*(copied_start + copied_end)

        # A pair can only reach the beam diameter if its midpoint distance is
        # no larger than the diameter plus both segment half-lengths.
        midpoint_difference = (
            current_midpoint[:, None, None, :]
            - copied_midpoint[None, :, :, :])
        midpoint_distance = np.linalg.norm(midpoint_difference, axis=3)
        candidate = midpoint_distance <= DIAMETER + maximum_length
        candidate[:, zero_shift, :] &= ~np.eye(
            number_of_elements, dtype=bool)
        first, shift_index, second = np.nonzero(candidate)

        # Exclude every element pair whose parent struts meet at the same
        # idealized lattice junction, including junctions across a periodic
        # cell boundary. Their finite-radius overlap belongs to the modeled
        # joint and must not be interpreted as beam-to-beam contact.
        first_start = strut_start[element_strut[first]]
        first_end = strut_end[element_strut[first]]
        second_start = (strut_start[element_strut[second]]
                        + integer_shifts[shift_index])
        second_end = (strut_end[element_strut[second]]
                      + integer_shifts[shift_index])
        junction_distance = np.minimum.reduce((
            np.linalg.norm(first_start - second_start, axis=1),
            np.linalg.norm(first_start - second_end, axis=1),
            np.linalg.norm(first_end - second_start, axis=1),
            np.linalg.norm(first_end - second_end, axis=1)))
        independent = junction_distance > 1.0e-8
        first = first[independent]
        shift_index = shift_index[independent]
        second = second[independent]

        current_distance = segment_distances(
            current_start[first], current_end[first],
            copied_start[shift_index, second],
            copied_end[shift_index, second])
        reference_shift = EDGE_LENGTH*integer_shifts[shift_index]
        reference_distance = segment_distances(
            reference_start[first], reference_end[first],
            reference_start[second] + reference_shift,
            reference_end[second] + reference_shift)

        if len(current_distance) == 0:
            minimum = np.inf
            local = None
        else:
            local = np.argmin(current_distance)
            minimum = current_distance[local]
        if local is not None and minimum < global_minimum:
            global_minimum = minimum
            closest_pair = {
                "frame": frame_index,
                "time": time,
                "first": int(first[local]),
                "second": int(second[local]),
                "shift": integer_shifts[shift_index[local]].copy(),
                "reference_distance": reference_distance[local],
                "current_distance": minimum,
                "reference_midpoint_first": 0.5*(
                    reference_start[first[local]]
                    + reference_end[first[local]]),
                "reference_midpoint_second": 0.5*(
                    reference_start[second[local]]
                    + reference_end[second[local]])
                    + reference_shift[local],
                "current_midpoint_first": current_midpoint[first[local]],
                "current_midpoint_second": copied_midpoint[
                    shift_index[local], second[local]],
            }
        measure = 0.4*time/STEP_DURATION
        output.append([time, measure, minimum])

        if ((frame_index + 1) % 50 == 0
                or frame_index + 1 == len(records)):
            print("{}: checked frame {} of {}".format(
                case_name, frame_index + 1, len(records)))

    return np.asarray(output), closest_pair


for case in ("compression", "tension", "shear_12"):
    result, closest_pair = minimum_clearances(case)
    path = os.path.join(
        HERE, "final_results", case, "contact_clearance.csv")
    np.savetxt(path, result, delimiter=",")
    print("\n{}".format(case))
    frame = np.argmin(result[:, 2])
    distance = result[frame, 2]
    print("  minimum independent-strut distance = {:.8e}".format(distance))
    print("  minimum beam-surface clearance      = {:.8e}".format(
        distance - DIAMETER))
    print("  deformation at minimum              = {:.6f}".format(
        result[frame, 1]))
    if closest_pair is None:
        print("  no independent strut pair entered the near-contact search "
              "radius")
    else:
        print("  elements = ({}, {}), image shift = {}, initial distance = "
              "{:.8e}".format(
                  closest_pair["first"], closest_pair["second"],
                  closest_pair["shift"],
                  closest_pair["reference_distance"]))
