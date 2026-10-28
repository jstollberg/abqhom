"""Geometry of the matrix-filled shell--beam sandwich unit cell."""

import gmsh


def _add_physical_group(dimension, entities, name):
    """Add and name one physical group."""
    group = gmsh.model.add_physical_group(dimension, entities)
    gmsh.model.set_physical_name(dimension, group, name)


def create_geometry(model_name, edge_length=1.0,
                    elements_per_strut=16):
    """Create the continuum, shell, and beam geometries.

    The continuum matrix occupies a cube centered at the origin. Its two
    faces normal to the second coordinate direction also carry shell
    elements. Eight beam struts connect the cell center to the vertices.
    Their corner nodes are shared with the continuum and shell meshes,
    whereas their interior nodes remain nonconforming and are subsequently
    coupled to the matrix by an embedded constraint.
    """
    if not gmsh.is_initialized():
        gmsh.initialize()
    if model_name in gmsh.model.list():
        gmsh.model.set_current(model_name)
        gmsh.model.remove()
    gmsh.model.add(model_name)

    gmsh.option.set_number("Geometry.Tolerance", 1.0e-8*edge_length)
    gmsh.option.set_number("Geometry.AutoCoherence", 0)

    half_edge = 0.5*edge_length
    geometry = gmsh.model.geo

    # Continuum cube. The ordering follows the first paper example. The
    # resulting surfaces are y-, y+, x-, z+, x+, and z-, respectively.
    coordinates = [(-half_edge, -half_edge, -half_edge),
                   (-half_edge, -half_edge,  half_edge),
                   ( half_edge, -half_edge,  half_edge),
                   ( half_edge, -half_edge, -half_edge),
                   (-half_edge,  half_edge, -half_edge),
                   (-half_edge,  half_edge,  half_edge),
                   ( half_edge,  half_edge,  half_edge),
                   ( half_edge,  half_edge, -half_edge)]
    matrix_points = [geometry.add_point(*point)
                     for point in coordinates]

    p1, p2, p3, p4, p5, p6, p7, p8 = matrix_points
    edges = [geometry.add_line(p1, p2),
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
             geometry.add_line(p4, p8)]

    surface_connectivity = [[edges[0], edges[1], edges[2], edges[3]],
                            [edges[4], edges[5], edges[6], edges[7]],
                            [edges[0], edges[9], edges[4], edges[8]],
                            [edges[1], edges[10], edges[5], edges[9]],
                            [edges[2], edges[11], edges[6], edges[10]],
                            [edges[3], edges[8], edges[7], edges[11]]]
    surfaces = []
    for connectivity in surface_connectivity:
        loop = geometry.add_curve_loop(connectivity, reorient=True)
        surfaces.append(geometry.add_plane_surface([loop]))

    surface_loop = geometry.add_surface_loop(surfaces)
    matrix_volume = geometry.add_volume([surface_loop])

    # BCC-type beam core. Reusing the matrix corner points makes the beam
    # endpoints conforming while retaining an independently meshed interior.
    beam_center = geometry.add_point(0.0, 0.0, 0.0)
    struts = [geometry.add_line(beam_center, vertex)
              for vertex in matrix_points]
    for strut in struts:
        geometry.mesh.set_transfinite_curve(
            strut, elements_per_strut + 1)

    geometry.synchronize()

    _add_physical_group(3, [matrix_volume], "CONTINUUM_MATRIX")
    _add_physical_group(2, [surfaces[0], surfaces[1]], "SHELL_FACES")
    _add_physical_group(1, struts, "BEAM_STRUTS")


def nodes_in_physical_group(dimension, name):
    """Return all nodes of a physical group, including its boundary."""
    groups = {
        gmsh.model.get_physical_name(dim, tag): (dim, tag)
        for dim, tag in gmsh.model.get_physical_groups()
    }
    if name not in groups:
        raise KeyError("Physical group '{}' does not exist.".format(name))
    dim, tag = groups[name]
    if dim != dimension:
        raise ValueError(
            "Physical group '{}' has dimension {}, expected {}."
            .format(name, dim, dimension))

    nodes = set()
    for entity in gmsh.model.get_entities_for_physical_group(dim, tag):
        entity_nodes = gmsh.model.mesh.get_nodes(
            dim, entity, includeBoundary=True)[0]
        nodes.update(int(node) for node in entity_nodes)
    return nodes

