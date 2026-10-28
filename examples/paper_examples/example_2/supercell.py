"""Geometry utilities for periodic Octet supercells."""

import os

import gmsh
import numpy as np

from abqhom.constants import PERIODIC_PAIRS


def _coordinate_key(coordinates):
    """Return a stable key for coincident lattice points."""
    return tuple(round(float(value), 12) for value in coordinates)


def create_geometry(model_name, edge_length, elements_per_strut,
                    cells_per_axis=2):
    """Create one conforming cubic Octet supercell.

    Points and struts shared by adjacent primitive cells are created only
    once. Consequently, internal cell interfaces are ordinary parts of the
    connected lattice and periodicity is imposed only on the outer boundary
    of the supercell.
    """
    if cells_per_axis < 1:
        raise ValueError("cells_per_axis must be a positive integer.")
    if not gmsh.is_initialized():
        gmsh.initialize()
    if model_name in gmsh.model.list():
        gmsh.model.set_current(model_name)
        gmsh.model.remove()
    gmsh.model.add(model_name)

    gmsh.option.set_number("Geometry.Tolerance", 1.0e-8*edge_length)
    gmsh.option.set_number("Mesh.ElementOrder", 1)
    geometry = gmsh.model.geo
    origin = -0.5*cells_per_axis*edge_length

    points = {}
    strut_keys = set()
    struts = []

    def point(coordinates):
        key = _coordinate_key(coordinates)
        if key not in points:
            points[key] = geometry.add_point(*coordinates)
        return key, points[key]

    def strut(first_coordinates, second_coordinates):
        first_key, first = point(first_coordinates)
        second_key, second = point(second_coordinates)
        key = tuple(sorted((first_key, second_key)))
        if key in strut_keys:
            return
        strut_keys.add(key)
        line = geometry.add_line(first, second)
        geometry.mesh.set_transfinite_curve(
            line, elements_per_strut + 1)
        struts.append(line)

    for cell_x in range(cells_per_axis):
        for cell_y in range(cells_per_axis):
            for cell_z in range(cells_per_axis):
                center = (
                    origin + (cell_x + 0.5)*edge_length,
                    origin + (cell_y + 0.5)*edge_length,
                    origin + (cell_z + 0.5)*edge_length)
                half_edge = 0.5*edge_length

                corners = {}
                for sign_x in (-1, 1):
                    for sign_y in (-1, 1):
                        for sign_z in (-1, 1):
                            signs = (sign_x, sign_y, sign_z)
                            corners[signs] = tuple(
                                center[axis] + half_edge*signs[axis]
                                for axis in range(3))

                face_centers = {}
                for axis in range(3):
                    for sign in (-1, 1):
                        coordinates = list(center)
                        coordinates[axis] += sign*half_edge
                        face_centers[(axis, sign)] = tuple(coordinates)

                for (axis, sign), face_center in face_centers.items():
                    for signs, corner in corners.items():
                        if signs[axis] == sign:
                            strut(face_center, corner)

                face_items = sorted(face_centers.items())
                for first_index, ((first_axis, _), first) in enumerate(
                        face_items):
                    for (second_axis, _), second in face_items[
                            first_index + 1:]:
                        if first_axis != second_axis:
                            strut(first, second)

    geometry.synchronize()
    group = gmsh.model.add_physical_group(1, struts)
    gmsh.model.set_physical_name(1, group, "OCTET_BEAMS")
    return len(points), len(struts)


def apply_eigenmode_imperfection(supercell, result_directory, mode_number,
                                 amplitude):
    """Perturb a supercell mesh by one normalized periodic eigenmode."""
    path = os.path.join(
        result_directory, "mode_{:02d}.csv".format(mode_number))
    data = np.loadtxt(path, delimiter=",")
    direction = {
        int(row[0]): np.asarray(row[1:4], dtype=float)
        for row in np.atleast_2d(data)}
    node_tags = [int(tag) for tag in supercell.node_sets["ALL_NODES"]]
    missing = [tag for tag in node_tags if tag not in direction]
    if missing:
        raise RuntimeError(
            "The buckling mode does not contain all Gmsh node tags.")

    normalization = max(np.linalg.norm(direction[tag])
                        for tag in node_tags)
    if normalization <= 0.0:
        raise RuntimeError("The selected buckling mode is zero.")
    direction = {tag: direction[tag]/normalization for tag in node_tags}

    residual = 0.0
    for slave, master in PERIODIC_PAIRS.items():
        if (slave not in supercell.node_sets
                or master not in supercell.node_sets):
            continue
        for slave_tag, master_tag in zip(supercell.node_sets[slave],
                                         supercell.node_sets[master]):
            mismatch = np.linalg.norm(
                direction[int(slave_tag)] - direction[int(master_tag)])
            residual = max(residual, mismatch)
    if residual > 1.0e-8:
        raise RuntimeError(
            "The imperfection is not periodic; residual = {:.3e}.".format(
                residual))

    for tag in node_tags:
        coordinates, parametric, _, _ = gmsh.model.mesh.get_node(tag)
        gmsh.model.mesh.set_node(
            tag, np.asarray(coordinates, dtype=float)
            + amplitude*direction[tag], parametric)
    return residual
