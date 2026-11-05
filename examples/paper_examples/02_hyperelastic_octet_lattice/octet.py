# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-or-later

"""Define the Octet topology, Gmsh geometry, and mode imperfection."""

import os

import numpy as np

from abqhom.utils import read_csv_file


def _coordinate_key(coordinates):
    """Return a stable key for coincident lattice points."""
    return tuple(round(float(value), 12) for value in coordinates)


def octet_topology(edge_length):
    """Return point coordinates and strut connectivity of the Octet cell."""
    half_edge = 0.5 * edge_length
    points = []

    # The primitive topology contains eight vertices and six face centers.
    corners = []
    for sign_x in (-1, 1):
        for sign_y in (-1, 1):
            for sign_z in (-1, 1):
                signs = (sign_x, sign_y, sign_z)
                corners.append((signs, len(points)))
                points.append(np.asarray(signs, dtype=float) * half_edge)

    face_centers = []
    for axis in range(3):
        for sign in (-1, 1):
            coordinates = np.zeros(3)
            coordinates[axis] = sign * half_edge
            face_centers.append(((axis, sign), len(points)))
            points.append(coordinates)

    # The 24 face-center-to-corner members and 12 face-center members form
    # the two periodic strut families of the conventional Octet cell.
    connectivity = []
    for (axis, sign), face_center in face_centers:
        for signs, corner in corners:
            if signs[axis] == sign:
                connectivity.append((face_center, corner))
    for first_index, ((first_axis, _), first) in enumerate(face_centers):
        for (second_axis, _), second in face_centers[first_index + 1:]:
            if first_axis != second_axis:
                connectivity.append((first, second))

    return (np.asarray(points, dtype=float),
            np.asarray(connectivity, dtype=int))


def supercell_topology(edge_length, cells_per_axis):
    """Return unique points and struts of a conforming Octet supercell."""
    if cells_per_axis < 1 or int(cells_per_axis) != cells_per_axis:
        raise ValueError("cells_per_axis must be a positive integer.")
    cells_per_axis = int(cells_per_axis)

    primitive_points, primitive_struts = octet_topology(edge_length)
    points = []
    point_indices = {}
    struts = []
    strut_keys = set()

    for cell_x in range(cells_per_axis):
        for cell_y in range(cells_per_axis):
            for cell_z in range(cells_per_axis):
                cell_indices = np.asarray(
                    [cell_x, cell_y, cell_z], dtype=float)
                center_index = 0.5 * (cells_per_axis - 1)
                offset = (cell_indices - center_index) * edge_length

                local_to_global = []
                for coordinates in primitive_points + offset:
                    # Shared points on internal cell interfaces receive one
                    # global index in the conforming supercell topology.
                    key = _coordinate_key(coordinates)
                    if key not in point_indices:
                        point_indices[key] = len(points)
                        points.append(coordinates)
                    local_to_global.append(point_indices[key])

                for first, second in primitive_struts:
                    pair = (local_to_global[first], local_to_global[second])
                    key = tuple(sorted(pair))
                    # The orientation-independent key likewise prevents
                    # duplicate struts on internal interfaces.
                    if key not in strut_keys:
                        strut_keys.add(key)
                        struts.append(pair)

    return (np.asarray(points, dtype=float),
            np.asarray(struts, dtype=int))


def create_geometry(model_name, edge_length, elements_per_strut,
                    cells_per_axis=1):
    """Create a conforming cubic Octet cell or supercell in Gmsh."""
    import gmsh

    if not gmsh.is_initialized():
        gmsh.initialize()
    if model_name in gmsh.model.list():
        gmsh.model.set_current(model_name)
        gmsh.model.remove()
    gmsh.model.add(model_name)

    gmsh.option.set_number("Geometry.Tolerance", 1.0e-8 * edge_length)
    gmsh.option.set_number("Mesh.ElementOrder", 1)

    geometry = gmsh.model.geo
    points, connectivity = supercell_topology(
        edge_length, cells_per_axis)
    point_tags = [geometry.add_point(*coordinates)
                  for coordinates in points]
    struts = [geometry.add_line(point_tags[first], point_tags[second])
              for first, second in connectivity]

    for strut in struts:
        geometry.mesh.set_transfinite_curve(
            strut, elements_per_strut + 1)
    geometry.synchronize()

    group = gmsh.model.add_physical_group(1, struts)
    gmsh.model.set_physical_name(1, group, "OCTET_BEAMS")


def _read_mode(path):
    """Read reference coordinates and translations of one buckling mode."""
    if not os.path.isfile(path):
        raise IOError(
            "The buckling-mode file '{}' does not exist. Run buckling.py "
            "before the finite-strain analyses.".format(path))

    data, _ = read_csv_file(path, has_header=True)
    data = np.atleast_2d(data)
    if data.shape[1] != 7 or not np.all(np.isfinite(data)):
        raise ValueError(
            "A mode file must contain the node tag, three reference "
            "coordinates, and three translational components. Rerun "
            "buckling.py to regenerate files in the current format.")

    tags = data[:, 0].astype(int)
    if not np.allclose(data[:, 0], tags) or len(set(tags)) != len(tags):
        raise ValueError("Mode-file node tags must be unique integers.")

    coordinates = {
        int(tag): np.asarray(values, dtype=float)
        for tag, values in zip(tags, data[:, 1:4])}
    translations = {
        int(tag): np.asarray(values, dtype=float)
        for tag, values in zip(tags, data[:, 4:7])}
    return coordinates, translations


def _repeat_primitive_mode(reference_coordinates, mode_translations,
                           current_nodes, edge_length, tolerance):
    """Map a primitive-cell mode periodically onto a compatible mesh."""
    if edge_length <= 0.0:
        raise ValueError("The primitive-cell edge length must be positive.")

    reference_values = np.asarray(list(reference_coordinates.values()))
    reference_minimum = np.min(reference_values, axis=0)
    reference_extent = np.ptp(reference_values, axis=0)
    if not np.allclose(
            reference_extent, edge_length,
            rtol=0.0, atol=1.0e-8 * edge_length):
        raise RuntimeError(
            "The stored mode does not span one cubic primitive cell.")

    current_values = np.asarray(
        [coordinates for coordinates, _ in current_nodes.values()])
    current_minimum = np.min(current_values, axis=0)
    cell_counts = np.ptp(current_values, axis=0) / edge_length
    if (np.any(cell_counts < 1.0 - 1.0e-8)
            or not np.allclose(
                cell_counts, np.round(cell_counts),
                rtol=0.0, atol=1.0e-8)):
        raise RuntimeError(
            "The target mesh does not span an integer number of primitive "
            "cells in every coordinate direction.")

    coordinate_tolerance = 1.0e-8 * edge_length

    def coordinate_key(coordinates):
        """Return a key at the coordinate-matching tolerance."""
        return tuple(int(round(float(value) / coordinate_tolerance))
                     for value in coordinates)

    def periodic_key(coordinates, minimum):
        """Return a coordinate key within the primitive reference cell."""
        phase = np.mod(coordinates - minimum, edge_length)
        near_boundary = np.minimum(
            np.abs(phase), np.abs(phase - edge_length))
        phase[near_boundary <= coordinate_tolerance] = 0.0
        return coordinate_key(reference_minimum + phase)

    periodic_mode = {}
    for tag, coordinates in reference_coordinates.items():
        key = periodic_key(coordinates, reference_minimum)
        if (key in periodic_mode
                and np.linalg.norm(
                    periodic_mode[key] - mode_translations[tag])
                > tolerance):
            raise RuntimeError(
                "The stored primitive-cell mode is not periodic.")
        periodic_mode[key] = mode_translations[tag]

    repeated_mode = {}
    for tag, (coordinates, _) in current_nodes.items():
        key = periodic_key(coordinates, current_minimum)
        if key not in periodic_mode:
            raise RuntimeError(
                "The target mesh is incompatible with the primitive-cell "
                "mode at node {}.".format(tag))
        repeated_mode[tag] = periodic_mode[key]
    return repeated_mode


def apply_eigenmode_imperfection(gmsh_cell, result_directory, amplitude,
                                 mode_number=1, tolerance=1.0e-10,
                                 primitive_edge_length=None):
    """Perturb a periodic mesh along a mode or its periodic repetition."""
    import gmsh

    from abqhom.constants import PERIODIC_PAIRS

    reference_coordinates, mode_translations = _read_mode(os.path.join(
        result_directory, "mode_{:02d}.csv".format(mode_number)))
    node_tags = [int(tag) for tag in gmsh_cell.node_sets["ALL_NODES"]]

    current_nodes = {}
    gmsh_cell.activate()
    for tag in node_tags:
        coordinates, parametric, _, _ = gmsh.model.mesh.get_node(tag)
        coordinates = np.asarray(coordinates, dtype=float)
        current_nodes[tag] = (coordinates, parametric)

    if primitive_edge_length is None:
        node_tag_set = set(node_tags)
        missing = [tag for tag in node_tags if tag not in mode_translations]
        extra = [tag for tag in mode_translations
                 if tag not in node_tag_set]
        if missing or extra:
            raise RuntimeError(
                "The buckling mode and current mesh have different node "
                "tags.")

        # Validate labels and reference coordinates before applying a mode
        # from the independently regenerated buckling model.
        coordinate_error = 0.0
        coordinate_scale = 1.0
        for tag in node_tags:
            coordinates = current_nodes[tag][0]
            coordinate_error = max(
                coordinate_error,
                np.linalg.norm(coordinates - reference_coordinates[tag]))
            coordinate_scale = max(
                coordinate_scale,
                np.linalg.norm(reference_coordinates[tag]))
        if coordinate_error > 1.0e-8 * coordinate_scale:
            raise RuntimeError(
                "The buckling mode was generated for a different mesh; "
                "maximum coordinate mismatch = {:.3e}.".format(
                    coordinate_error))
    else:
        mode_translations = _repeat_primitive_mode(
            reference_coordinates, mode_translations, current_nodes,
            primitive_edge_length, tolerance)

    normalization = max(
        np.linalg.norm(mode_translations[tag]) for tag in node_tags)
    if normalization <= 0.0:
        raise RuntimeError("The selected buckling mode is zero.")

    # The prescribed amplitude is therefore the largest nodal perturbation.
    direction = {
        tag: mode_translations[tag] / normalization for tag in node_tags}

    # The eigenproblem must produce matching translations on periodic pairs
    residual = 0.0
    for secondary, primary in PERIODIC_PAIRS.items():
        sets = gmsh_cell.node_sets
        pair_available = secondary in sets and primary in sets
        if not pair_available:
            continue
        secondary_tags = sets[secondary]
        primary_tags = sets[primary]
        if len(secondary_tags) != len(primary_tags):
            raise RuntimeError(
                "Periodic node sets do not form one-to-one pairs.")
        for secondary_tag, primary_tag in zip(secondary_tags, primary_tags):
            mismatch = np.linalg.norm(
                direction[int(secondary_tag)] - direction[int(primary_tag)])
            residual = max(residual, mismatch)

    if residual > tolerance:
        raise RuntimeError(
            "The imperfection is not periodic; residual = {:.3e}.".format(
                residual))

    # Only translational mode components alter the Gmsh nodal positions;
    # finite beam rotations remain unknowns of the subsequent analysis.
    for tag in node_tags:
        coordinates, parametric = current_nodes[tag]
        coordinates = coordinates.copy()
        coordinates += amplitude * direction[tag]
        gmsh.model.mesh.set_node(tag, coordinates, parametric)

    return residual
