"""
Geometry fixtures shared by the abqhom unit tests.

Each builder creates a small Gmsh model with one or more physical groups
and returns the model name. They only use Gmsh, so they can be imported
without Abaqus.

"""
import gmsh


def _new_model(name):
    if not gmsh.is_initialized():
        gmsh.initialize()
    if name in gmsh.model.list():
        gmsh.model.set_current(name)
        gmsh.model.remove()
    gmsh.model.add(name)
    return gmsh.model.geo


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
    """The four edge struts of a square (2d); one group 'STRUTS'."""
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


def plus_lattice_2d(name="PLUS", size=1.0):
    """A 2d cell whose struts cross the boundary at the edge midpoints,
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
