"""Geometry and periodic eigenmode utilities for the Octet example."""

import os

import gmsh
import numpy as np

from abqhom.constants import PERIODIC_PAIRS


def create_geometry(model_name, edge_length, elements_per_strut):
    """Create a conventional cubic representation of the Octet lattice."""
    if not gmsh.is_initialized():
        gmsh.initialize()
    if model_name in gmsh.model.list():
        gmsh.model.set_current(model_name)
        gmsh.model.remove()
    gmsh.model.add(model_name)

    gmsh.option.set_number("Geometry.Tolerance", 1.0e-8*edge_length)
    gmsh.option.set_number("Mesh.ElementOrder", 1)
    geometry = gmsh.model.geo
    half_edge = 0.5*edge_length

    corners = {}
    for sign_x in (-1, 1):
        for sign_y in (-1, 1):
            for sign_z in (-1, 1):
                signs = (sign_x, sign_y, sign_z)
                corners[signs] = geometry.add_point(
                    sign_x*half_edge, sign_y*half_edge,
                    sign_z*half_edge)

    face_centers = {}
    for axis in range(3):
        for sign in (-1, 1):
            coordinates = [0.0, 0.0, 0.0]
            coordinates[axis] = sign*half_edge
            face_centers[(axis, sign)] = geometry.add_point(*coordinates)

    struts = []
    for (axis, sign), face_center in face_centers.items():
        for signs, corner in corners.items():
            if signs[axis] == sign:
                struts.append(geometry.add_line(face_center, corner))

    face_items = sorted(face_centers.items())
    for first_index, ((first_axis, _), first) in enumerate(face_items):
        for (second_axis, _), second in face_items[first_index + 1:]:
            if first_axis != second_axis:
                struts.append(geometry.add_line(first, second))

    for strut in struts:
        geometry.mesh.set_transfinite_curve(
            strut, elements_per_strut + 1)
    geometry.synchronize()

    group = gmsh.model.add_physical_group(1, struts)
    gmsh.model.set_physical_name(1, group, "OCTET_BEAMS")


def read_mode(path):
    """Read a normalized translational buckling mode by node tag."""
    data = np.loadtxt(path, delimiter=",")
    return {int(row[0]): np.asarray(row[1:4], dtype=float)
            for row in np.atleast_2d(data)}


def apply_eigenmode_imperfection(unit_cell, result_directory, amplitude,
                                 angle=0.0):
    """Perturb a periodic mesh in its lowest, two-dimensional eigenspace."""
    mode_1 = read_mode(os.path.join(result_directory, "mode_01.csv"))
    mode_2 = read_mode(os.path.join(result_directory, "mode_02.csv"))
    node_tags = [int(tag) for tag in unit_cell.node_sets["ALL_NODES"]]

    missing = [tag for tag in node_tags
               if tag not in mode_1 or tag not in mode_2]
    if missing:
        raise RuntimeError(
            "The buckling modes do not contain all Gmsh node tags.")

    direction = {
        tag: (np.cos(angle)*mode_1[tag]
              + np.sin(angle)*mode_2[tag])
        for tag in node_tags}
    normalization = max(np.linalg.norm(value)
                        for value in direction.values())
    if normalization <= 0.0:
        raise RuntimeError("The selected buckling-mode combination is zero.")
    direction = {tag: value/normalization
                 for tag, value in direction.items()}

    residual = 0.0
    for slave, master in PERIODIC_PAIRS.items():
        if (slave not in unit_cell.node_sets
                or master not in unit_cell.node_sets):
            continue
        for slave_tag, master_tag in zip(unit_cell.node_sets[slave],
                                         unit_cell.node_sets[master]):
            mismatch = np.linalg.norm(
                direction[int(slave_tag)] - direction[int(master_tag)])
            residual = max(residual, mismatch)
    if residual > 1.0e-10:
        raise RuntimeError(
            "The imperfection is not periodic; residual = {:.3e}.".format(
                residual))

    for tag in node_tags:
        coordinates, parametric, _, _ = gmsh.model.mesh.get_node(tag)
        gmsh.model.mesh.set_node(
            tag, np.asarray(coordinates, dtype=float)
            + amplitude*direction[tag], parametric)

    return residual
