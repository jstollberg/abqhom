# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-or-later

"""Define the Gmsh geometry of the elastoplastic sandwich cell."""


def create_geometry(model_name, edge_length, elements_per_strut):
    """Create the continuum matrix, shell skins, and BCC beam core."""
    import gmsh

    def add_physical_group(dimension, entities, name):
        """Add and name one physical group."""
        group = gmsh.model.add_physical_group(dimension, entities)
        gmsh.model.set_physical_name(dimension, group, name)

    if not gmsh.is_initialized():
        gmsh.initialize()
    if model_name in gmsh.model.list():
        gmsh.model.set_current(model_name)
        gmsh.model.remove()
    gmsh.model.add(model_name)

    gmsh.option.set_number("Geometry.Tolerance", 1.0e-8 * edge_length)

    # Prevent automatic merging of the independently meshed beam interiors
    # with the continuum while retaining the explicitly shared endpoints.
    gmsh.option.set_number("Geometry.AutoCoherence", 0)

    half_edge = 0.5 * edge_length
    geometry = gmsh.model.geo

    # Four corners at y=-h, followed by the corresponding corners at y=+h.
    corner_coordinates = [
        (-half_edge, -half_edge, -half_edge),
        (-half_edge, -half_edge, half_edge),
        (half_edge, -half_edge, half_edge),
        (half_edge, -half_edge, -half_edge),
        (-half_edge, half_edge, -half_edge),
        (-half_edge, half_edge, half_edge),
        (half_edge, half_edge, half_edge),
        (half_edge, half_edge, -half_edge),
    ]
    points = [geometry.add_point(*coordinates)
              for coordinates in corner_coordinates]
    p1, p2, p3, p4, p5, p6, p7, p8 = points

    edges = [
        geometry.add_line(p1, p2),
        geometry.add_line(p2, p3),
        geometry.add_line(p3, p4),
        geometry.add_line(p4, p1),
        geometry.add_line(p5, p6),
        geometry.add_line(p6, p7),
        geometry.add_line(p7, p8),
        geometry.add_line(p8, p5),
        geometry.add_line(p1, p5),
        geometry.add_line(p2, p6),
        geometry.add_line(p3, p7),
        geometry.add_line(p4, p8),
    ]

    # Surface order: y-, y+, x-, z+, x+, z-.
    surface_connectivity = [
        [edges[0], edges[1], edges[2], edges[3]],
        [edges[4], edges[5], edges[6], edges[7]],
        [edges[0], edges[9], edges[4], edges[8]],
        [edges[1], edges[10], edges[5], edges[9]],
        [edges[2], edges[11], edges[6], edges[10]],
        [edges[3], edges[8], edges[7], edges[11]],
    ]
    surfaces = []
    for connectivity in surface_connectivity:
        loop = geometry.add_curve_loop(connectivity, reorient=True)
        surfaces.append(geometry.add_plane_surface([loop]))

    surface_loop = geometry.add_surface_loop(surfaces)
    matrix = geometry.add_volume([surface_loop])

    # Reusing the matrix vertices makes the beam endpoints conforming. The
    # remaining beam nodes are subsequently coupled by an embedded region.
    center = geometry.add_point(0.0, 0.0, 0.0)
    struts = [geometry.add_line(center, vertex) for vertex in points]
    for strut in struts:
        geometry.mesh.set_transfinite_curve(
            strut, elements_per_strut + 1)

    geometry.synchronize()
    add_physical_group(3, [matrix], "CONTINUUM_MATRIX")
    add_physical_group(
        2, [surfaces[0], surfaces[1]], "SHELL_FACES")
    add_physical_group(1, struts, "BEAM_STRUTS")
