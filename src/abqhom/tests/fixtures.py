"""
Geometry fixtures shared by the abqhom unit tests.

Each builder creates a small Gmsh model with one or more physical groups
and returns the model name. They only use Gmsh, so they can be imported
without Abaqus.

"""
import gmsh


def _new_model(name):
    """Create an empty Gmsh model, replacing one with the same name."""
    if not gmsh.is_initialized():
        gmsh.initialize()
    if name in gmsh.model.list():
        gmsh.model.set_current(name)
        gmsh.model.remove()
    gmsh.model.add(name)
    return gmsh.model.geo


def discrete_periodic_plates(name="PLATES", size=1.0, perturb=0.0):
    """Two parallel plates at z=0 and z=size, each a 3x3 grid of nodes, as
    a pre-built discrete surface mesh with one physical group 'PLATES'.
    Used to exercise GmshRVE.read_from_mesh. With perturb=0 the mesh is
    node-matched (periodic); a non-zero perturb shifts the top plate's
    centre node in x so that the mesh is no longer periodic."""
    if not gmsh.is_initialized():
        gmsh.initialize()
    if name in gmsh.model.list():
        gmsh.model.set_current(name)
        gmsh.model.remove()
    gmsh.model.add(name)

    grid = [0.0, 0.5*size, size]
    nodes = []
    for z in (0.0, size):
        for x in grid:
            for y in grid:
                nodes.append((x, y, z))
    if perturb != 0.0:
        # Move the top plate's centre node (index 9 + 1*3 + 1) off the grid.
        nodes[13] = (0.5*size + perturb, 0.5*size, size)

    tris = []
    for plate in (0, 1):
        base = plate*9
        for ix in range(2):
            for iy in range(2):
                a = base + ix*3 + iy
                b = base + (ix + 1)*3 + iy
                c = base + (ix + 1)*3 + (iy + 1)
                d = base + ix*3 + (iy + 1)
                tris.append((a, b, c))
                tris.append((a, c, d))

    surface = gmsh.model.add_discrete_entity(2)
    coords = [v for node in nodes for v in node]
    gmsh.model.mesh.add_nodes(2, surface, list(range(1, len(nodes) + 1)),
                              coords)
    conn = [t + 1 for tri in tris for t in tri]
    gmsh.model.mesh.add_elements_by_type(surface, 2,
                                         list(range(1, len(tris) + 1)), conn)
    group = gmsh.model.add_physical_group(2, [surface])
    gmsh.model.set_physical_name(2, group, "PLATES")
    return name


def discrete_periodic_square(name="SQUARE", size=1.0, perturb=0.0):
    """A pre-built 3x3 triangular mesh of a square in the x-y plane.

    With ``perturb=0`` the nodes on opposite edges are matched. A non-zero
    perturb shifts the midpoint of the right edge in y and makes the mesh
    non-periodic without changing its bounding box.
    """
    if not gmsh.is_initialized():
        gmsh.initialize()
    if name in gmsh.model.list():
        gmsh.model.set_current(name)
        gmsh.model.remove()
    gmsh.model.add(name)

    grid = [0.0, 0.5*size, size]
    nodes = [(x, y, 0.0) for x in grid for y in grid]
    if perturb != 0.0:
        nodes[7] = (size, 0.5*size + perturb, 0.0)

    tris = []
    for ix in range(2):
        for iy in range(2):
            a = ix*3 + iy
            b = (ix + 1)*3 + iy
            c = (ix + 1)*3 + (iy + 1)
            d = ix*3 + (iy + 1)
            tris.append((a, b, c))
            tris.append((a, c, d))

    surface = gmsh.model.add_discrete_entity(2)
    coords = [value for node in nodes for value in node]
    gmsh.model.mesh.add_nodes(2, surface, list(range(1, len(nodes) + 1)),
                              coords)
    connectivity = [node + 1 for tri in tris for node in tri]
    gmsh.model.mesh.add_elements_by_type(
        surface, 2, list(range(1, len(tris) + 1)), connectivity)
    group = gmsh.model.add_physical_group(2, [surface])
    gmsh.model.set_physical_name(2, group, "SOLID")
    return name


def simple_cubic(name="SIMPLE_CUBIC", length=1.0):
    """Twelve edge struts of a cube; one physical group 'STRUTS'."""
    g = _new_model(name)
    h = 0.5*length
    coords = [(x, y, z) for x in (-h, h) for y in (-h, h) for z in (-h, h)]
    corners = [g.add_point(c[0], c[1], c[2]) for c in coords]
    struts = []
    for i in range(8):
        for j in range(i + 1, 8):
            shared = 0
            for a, b in zip(coords[i], coords[j]):
                if a == b:
                    shared += 1
            if shared == 2:
                struts.append(g.add_line(corners[i], corners[j]))
    g.synchronize()
    group = gmsh.model.add_physical_group(1, struts)
    gmsh.model.set_physical_name(1, group, "STRUTS")
    return name


def bcc(name="BCC", length=1.0):
    """Cubic edge struts ('CUBIC') plus body diagonals ('BCC')."""
    g = _new_model(name)
    h = 0.5*length
    coords = [(x, y, z) for x in (-h, h) for y in (-h, h) for z in (-h, h)]
    corners = [g.add_point(c[0], c[1], c[2]) for c in coords]
    center = g.add_point(0.0, 0.0, 0.0)
    cubic = []
    for i in range(8):
        for j in range(i + 1, 8):
            shared = 0
            for a, b in zip(coords[i], coords[j]):
                if a == b:
                    shared += 1
            if shared == 2:
                cubic.append(g.add_line(corners[i], corners[j]))
    diagonals = [g.add_line(center, corner) for corner in corners]
    g.synchronize()
    group = gmsh.model.add_physical_group(1, cubic)
    gmsh.model.set_physical_name(1, group, "CUBIC")
    group = gmsh.model.add_physical_group(1, diagonals)
    gmsh.model.set_physical_name(1, group, "BCC")
    return name


def square_lattice(name="SQUARE", size=1.0):
    """The four edge struts of a square (2D); one group 'STRUTS'."""
    g = _new_model(name)
    h = 0.5*size
    p = [g.add_point(-h, -h, 0.0), g.add_point(h, -h, 0.0),
         g.add_point(h, h, 0.0), g.add_point(-h, h, 0.0)]
    struts = [g.add_line(p[0], p[1]), g.add_line(p[1], p[2]),
              g.add_line(p[2], p[3]), g.add_line(p[3], p[0])]
    g.synchronize()
    group = gmsh.model.add_physical_group(1, struts)
    gmsh.model.set_physical_name(1, group, "STRUTS")
    return name


def solid_cube(name="SOLID", length=1.0):
    """A solid continuum cube as one physical volume 'SOLID'."""
    if not gmsh.is_initialized():
        gmsh.initialize()
    if name in gmsh.model.list():
        gmsh.model.set_current(name)
        gmsh.model.remove()
    gmsh.model.add(name)
    h = 0.5*length
    box = gmsh.model.occ.add_box(-h, -h, -h, length, length, length)
    gmsh.model.occ.synchronize()
    group = gmsh.model.add_physical_group(3, [box])
    gmsh.model.set_physical_name(3, group, "SOLID")
    return name


def parallel_plates_3d(name="PLATES", size=1.0, spacing=1.0):
    """Two flat shell plates at z = +/- spacing/2 spanning a size x size
    square; one surface group 'SHELLS'. A periodic stack of membranes with
    a 3D bounding box, so the surface elements map to shells."""
    g = _new_model(name)
    a = 0.5*size
    hz = 0.5*spacing
    surfaces = []
    for z in (-hz, hz):
        pts = [g.add_point(-a, -a, z), g.add_point(a, -a, z),
               g.add_point(a, a, z), g.add_point(-a, a, z)]
        lines = [g.add_line(pts[0], pts[1]), g.add_line(pts[1], pts[2]),
                 g.add_line(pts[2], pts[3]), g.add_line(pts[3], pts[0])]
        loop = g.add_curve_loop(lines)
        surfaces.append(g.add_plane_surface([loop]))
    g.synchronize()
    group = gmsh.model.add_physical_group(2, surfaces)
    gmsh.model.set_physical_name(2, group, "SHELLS")
    return name


def mixed_dimensional_cube(name="MIXED", size=1.0, nodes_per_edge=3):
    """A conforming cube with continuum, shell, and beam constituents."""
    g = _new_model(name)
    h = 0.5*size
    coordinates = [(-h, -h, -h), (-h, -h, h),
                   (h, -h, h), (h, -h, -h),
                   (-h, h, -h), (-h, h, h),
                   (h, h, h), (h, h, -h)]
    points = [g.add_point(*coordinate) for coordinate in coordinates]
    line_connectivity = [(0, 1), (1, 2), (2, 3), (3, 0),
                         (4, 5), (5, 6), (6, 7), (7, 4),
                         (0, 4), (1, 5), (2, 6), (3, 7)]
    lines = [g.add_line(points[first], points[second])
             for first, second in line_connectivity]
    surface_connectivity = [(0, 1, 2, 3), (4, 5, 6, 7),
                            (0, 9, 4, 8), (1, 10, 5, 9),
                            (2, 11, 6, 10), (3, 8, 7, 11)]
    surfaces = []
    for connectivity in surface_connectivity:
        loop = g.add_curve_loop([lines[index] for index in connectivity],
                                reorient=True)
        surfaces.append(g.add_plane_surface([loop]))
    volume = g.add_volume([g.add_surface_loop(surfaces)])

    for line in lines:
        g.mesh.set_transfinite_curve(line, nodes_per_edge)
    for surface in surfaces:
        g.mesh.set_transfinite_surface(surface)
        g.mesh.set_recombine(2, surface)
    g.mesh.set_transfinite_volume(volume)
    g.synchronize()

    group = gmsh.model.add_physical_group(3, [volume])
    gmsh.model.set_physical_name(3, group, "CONTINUUM")
    group = gmsh.model.add_physical_group(2, surfaces)
    gmsh.model.set_physical_name(2, group, "SHELLS")
    group = gmsh.model.add_physical_group(1, lines)
    gmsh.model.set_physical_name(1, group, "BEAMS")
    return name


def matrix_with_boundary_fiber(name="MATRIXFIBER", size=1.0, n=5):
    """A continuum cube (transfinite hex, group 'MATRIX') with a single beam
    fiber running the full height, so its end nodes lie ON the z-faces
    (group 'FIBER'). Meant to be tied to the matrix with an embedded-region
    constraint: its boundary nodes exercise the exclusion of embedded
    (translation-eliminated) nodes from the translational periodic
    equations."""
    g = _new_model(name)
    h = 0.5*size
    c = [(-h, -h, -h), (-h, -h, h), (h, -h, h), (h, -h, -h),
         (-h, h, -h), (-h, h, h), (h, h, h), (h, h, -h)]
    pts = [g.add_point(*p) for p in c]
    lc = [[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6, 7], [7, 4],
          [0, 4], [1, 5], [2, 6], [3, 7]]
    lines = [g.add_line(pts[a], pts[b]) for a, b in lc]
    sc = [[0, 1, 2, 3], [4, 5, 6, 7], [0, 9, 4, 8], [1, 10, 5, 9],
          [2, 11, 6, 10], [3, 8, 7, 11]]
    surfaces = []
    for s in sc:
        loop = g.add_curve_loop([lines[i] for i in s], reorient=True)
        surfaces.append(g.add_plane_surface([loop]))
    vol = g.add_volume([g.add_surface_loop(surfaces)])
    for i in lines:
        g.mesh.set_transfinite_curve(i, n)
    for i in surfaces:
        g.mesh.set_transfinite_surface(i)
        g.mesh.set_recombine(2, i)
    g.mesh.set_transfinite_volume(vol)
    fp1 = g.add_point(-0.25*size, 0.25*size, -h)
    fp2 = g.add_point(-0.25*size, 0.25*size, h)
    fiber = g.add_line(fp1, fp2)
    g.mesh.set_transfinite_curve(fiber, 3)
    g.synchronize()
    gm = gmsh.model.add_physical_group(3, [vol])
    gmsh.model.set_physical_name(3, gm, "MATRIX")
    gf = gmsh.model.add_physical_group(1, [fiber])
    gmsh.model.set_physical_name(1, gf, "FIBER")
    return name


def plus_lattice_2d(name="PLUS", size=1.0):
    """A 2D cell whose struts cross the boundary at the edge midpoints,
    i.e. not at the corners; one group 'STRUTS'."""
    g = _new_model(name)
    h = 0.5*size
    center = g.add_point(0.0, 0.0, 0.0)
    mids = [g.add_point(0.0, -h, 0.0), g.add_point(h, 0.0, 0.0),
            g.add_point(0.0, h, 0.0), g.add_point(-h, 0.0, 0.0)]
    struts = [g.add_line(center, m) for m in mids]
    g.synchronize()
    group = gmsh.model.add_physical_group(1, struts)
    gmsh.model.set_physical_name(1, group, "STRUTS")
    return name
