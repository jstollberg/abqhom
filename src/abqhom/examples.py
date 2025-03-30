import gmsh
from abqhom.RVE import (add_boundary_groups, add_periodic_constraints,
                        create_boundary_sets)

def fancy_RVE_2d(model_name="RVE", x=0.0, y=0.0, dx=1.0, dy=1.0, radius=0.25,
                 lc=None, gui=False):
    # initialize gmsh and add a new model
    if not gmsh.isInitialized():
        gmsh.initialize()
    gmsh.model.add(model_name)

    # create geometry
    gmsh.model.occ.addRectangle(x - dx/2, y - dy/2, 0.0, dx, dy)
    gmsh.model.occ.addDisk(x - dx/2, y, 0.0, radius, radius)
    gmsh.model.occ.addDisk(x + dx/2, y, 0.0, radius, radius)
    gmsh.model.occ.addDisk(x, y - dy/2, 0.0, radius, radius)
    gmsh.model.occ.addDisk(x, y + dy/2, 0.0, radius, radius)
    gmsh.model.occ.cut([(2, 1)], [(2, 2)])
    gmsh.model.occ.cut([(2, 1)], [(2, 3)])
    gmsh.model.occ.cut([(2, 1)], [(2, 4)])
    gmsh.model.occ.cut([(2, 1)], [(2, 5)])
    gmsh.model.occ.synchronize()

    # create map containing information on mesh periodicity and physical groups
    group_map = {"POINT_1": [13],
                 "POINT_2": [10],
                 "POINT_3": [7],
                 "POINT_4": [4],
                 "EDGE_1": [12, 10],
                 "EDGE_2": [9, 7],
                 "EDGE_3": [4, 6],
                 "EDGE_4": [13, 3]}

    add_boundary_groups(model_name, group_map)
    add_periodic_constraints(model_name, group_map, dx, dy)

    # set mesh size
    if lc is not None:
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeMin", lc)
        gmsh.option.setNumber("Mesh.MeshSizeMax", lc)

    # create mesh using Delauny algorithm
    gmsh.option.setNumber("Mesh.Algorithm", 5)
    gmsh.model.mesh.generate(2)

    if gui:
        gmsh.fltk.run()

    return model_name, group_map

def fancy_RVE_3d(model_name="RVE", x=0.0, y=0.0, z=0.0, dx=1.0, dy=1.0, dz=1.0,
                 radius=0.25, lc=None, gui=False):
    # initialize gmsh and add a new model
    if not gmsh.isInitialized():
        gmsh.initialize()
    gmsh.model.add(model_name)

    # create geometry
    gmsh.model.occ.addBox(x - dx/2, y - dy/2, z - dz/2, dx, dy, dz)
    gmsh.model.occ.addCylinder(x - dx/2, y, z - dz/2, 0.0, 0.0, dz, radius)
    gmsh.model.occ.addCylinder(x + dx/2, y, z - dz/2, 0.0, 0.0, dz, radius)
    gmsh.model.occ.addCylinder(x, y - dy/2, z - dz/2, 0.0, 0.0, dz, radius)
    gmsh.model.occ.addCylinder(x, y + dy/2, z - dz/2, 0.0, 0.0, dz, radius)
    gmsh.model.occ.cut([(3, 1)], [(3, 2)])
    gmsh.model.occ.cut([(3, 1)], [(3, 3)])
    gmsh.model.occ.cut([(3, 1)], [(3, 4)])
    gmsh.model.occ.cut([(3, 1)], [(3, 5)])
    gmsh.model.occ.synchronize()

    # create map containing information on mesh periodicity and physical groups
    group_map = {"POINT_1": [1],
                 "POINT_2": [4],
                 "POINT_3": [18],
                 "POINT_4": [7],
                 "POINT_5": [13],
                 "POINT_6": [24],
                 "POINT_7": [21],
                 "POINT_8": [10],
                 "EDGE_1": [4],
                 "EDGE_2": [3, 20],
                 "EDGE_3": [32],
                 "EDGE_4": [1, 7],
                 "EDGE_5": [38],
                 "EDGE_6": [26, 24],
                 "EDGE_7": [35],
                 "EDGE_8": [13, 11],
                 "EDGE_9": [27, 29],
                 "EDGE_10": [14, 16],
                 "EDGE_11": [23, 21],
                 "EDGE_12": [10, 8],
                 "FACE_1": [1, 7],
                 "FACE_2": [13, 11],
                 "FACE_3": [5, 14],
                 "FACE_4": [8, 10],
                 "FACE_5": [4],
                 "FACE_6": [2]}

    add_boundary_groups(model_name, group_map)
    add_periodic_constraints(model_name, group_map, dx, dy, dz)

    # set mesh size
    if lc is not None:
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeMin", lc)
        gmsh.option.setNumber("Mesh.MeshSizeMax", lc)

    # create mesh using Delauny algorithm
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)
    gmsh.model.mesh.generate(3)

    if gui:
        gmsh.fltk.run()

    return model_name, group_map

def hole_RVE_2d(model_name="RVE", x=0.0, y=0.0, dx=1.0, dy=1.0, radius=0.25,
                lc=None, gui=False):
    # initialize gmsh and add a new model
    if not gmsh.isInitialized():
        gmsh.initialize()
    gmsh.model.add(model_name)

    # create geometry
    gmsh.model.occ.addRectangle(x - dx/2, y - dy/2, 0.0, dx, dy)
    gmsh.model.occ.addDisk(x, y, 0.0, radius, radius)
    gmsh.model.occ.cut([(2, 1)], [(2, 2)])
    gmsh.model.occ.synchronize()

    # create map containing information on mesh periodicity and physical groups
    group_map = {"POINT_1": [6],
                 "POINT_2": [7],
                 "POINT_3": [9],
                 "POINT_4": [8],
                 "EDGE_1": [6],
                 "EDGE_2": [8],
                 "EDGE_3": [9],
                 "EDGE_4": [7]}

    add_boundary_groups(model_name, group_map)
    add_periodic_constraints(model_name, group_map, dx, dy)

    # set mesh size
    if lc is not None:
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeMin", lc)
        gmsh.option.setNumber("Mesh.MeshSizeMax", lc)

    # create mesh using Delauny algorithm
    gmsh.option.setNumber("Mesh.Algorithm", 5)
    gmsh.model.mesh.generate(2)

    if gui:
        gmsh.fltk.run()

    return model_name, group_map

def hole_RVE_3d(model_name="RVE", x=0.0, y=0.0, z=0.0, dx=1.0, dy=1.0, dz=1.0,
                radius=0.25, lc=None, gui=False):
    # initialize gmsh and add a new model
    if not gmsh.isInitialized():
        gmsh.initialize()
    gmsh.model.add(model_name)

    # create geometry
    gmsh.model.occ.addBox(x - dx/2, y - dy/2, z - dz/2, dx, dy, dz)
    gmsh.model.occ.addCylinder(x, y, z - dz/2, 0.0, 0.0, dz, radius)
    gmsh.model.occ.cut([(3, 1)], [(3, 2)])
    gmsh.model.occ.synchronize()

    # create map containing information on mesh periodicity and physical groups
    group_map = {"POINT_1": [1],
                 "POINT_2": [2],
                 "POINT_3": [4],
                 "POINT_4": [3],
                 "POINT_5": [6],
                 "POINT_6": [5],
                 "POINT_7": [9],
                 "POINT_8": [7],
                 "EDGE_1": [1],
                 "EDGE_2": [4],
                 "EDGE_3": [3],
                 "EDGE_4": [2],
                 "EDGE_5": [6],
                 "EDGE_6": [13],
                 "EDGE_7": [12],
                 "EDGE_8": [9],
                 "EDGE_9": [5],
                 "EDGE_10": [7],
                 "EDGE_11": [11],
                 "EDGE_12": [8],
                 "FACE_1": [1],
                 "FACE_2": [6],
                 "FACE_3": [2],
                 "FACE_4": [4],
                 "FACE_5": [5],
                 "FACE_6": [3]}

    add_boundary_groups(model_name, group_map)
    add_periodic_constraints(model_name, group_map, dx, dy, dz)

    # set mesh size
    if lc is not None:
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeMin", lc)
        gmsh.option.setNumber("Mesh.MeshSizeMax", lc)

    # create mesh using Delauny algorithm
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)
    gmsh.model.mesh.generate(3)

    if gui:
        gmsh.fltk.run()

    return model_name, group_map

def simple_RVE_2d(model_name="RVE", x=0.0, y=0.0, dx=1.0, dy=1.0, lc=None,
                  gui=False):
    # initialize gmsh and add a new model
    if not gmsh.isInitialized():
        gmsh.initialize()
    gmsh.model.add(model_name)

    # create geometry
    gmsh.model.occ.addRectangle(x - dx/2, y - dy/2, 0.0, dx, dy)
    gmsh.model.occ.synchronize()

    # create map containing information on mesh periodicity and physical groups
    group_map = {"POINT_1": [1],
                 "POINT_2": [2],
                 "POINT_3": [3],
                 "POINT_4": [4],
                 "EDGE_1": [1],
                 "EDGE_2": [2],
                 "EDGE_3": [3],
                 "EDGE_4": [4]}

    add_boundary_groups(model_name, group_map)
    add_periodic_constraints(model_name, group_map, dx, dy)

    # set mesh size
    if lc is not None:
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeMin", lc)
        gmsh.option.setNumber("Mesh.MeshSizeMax", lc)

    # create mesh using Delauny algorithm
    gmsh.option.setNumber("Mesh.Algorithm", 5)
    gmsh.model.mesh.generate(2)

    if gui:
        gmsh.fltk.run()

    return model_name, group_map

def simple_RVE_3d(model_name="RVE", x=0.0, y=0.0, z=0.0, dx=1.0, dy=1.0,
                  dz=1.0, lc=None, gui=False):
    # initialize gmsh and add a new model
    if not gmsh.isInitialized():
        gmsh.initialize()
    gmsh.model.add(model_name)

    # create geometry
    gmsh.model.occ.addBox(x - dx/2, y - dy/2, z - dz/2, dx, dy, dz)
    gmsh.model.occ.synchronize()

    # create map containing information on mesh periodicity and physical groups
    group_map = {"POINT_1": [1],
                 "POINT_2": [2],
                 "POINT_3": [4],
                 "POINT_4": [3],
                 "POINT_5": [5],
                 "POINT_6": [6],
                 "POINT_7": [8],
                 "POINT_8": [7],
                 "EDGE_1": [1],
                 "EDGE_2": [4],
                 "EDGE_3": [3],
                 "EDGE_4": [2],
                 "EDGE_5": [5],
                 "EDGE_6": [8],
                 "EDGE_7": [7],
                 "EDGE_8": [6],
                 "EDGE_9": [9],
                 "EDGE_10": [10],
                 "EDGE_11": [11],
                 "EDGE_12": [12],
                 "FACE_1": [1],
                 "FACE_2": [2],
                 "FACE_3": [3],
                 "FACE_4": [4],
                 "FACE_5": [5],
                 "FACE_6": [6]}

    add_boundary_groups(model_name, group_map)
    add_periodic_constraints(model_name, group_map, dx, dy, dz)

    # set mesh size
    if lc is not None:
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeMin", lc)
        gmsh.option.setNumber("Mesh.MeshSizeMax", lc)

    # create mesh using Delauny algorithm
    gmsh.option.setNumber("Mesh.Algorithm", 5)
    gmsh.model.mesh.generate(3)

    if gui:
        gmsh.fltk.run()

    return model_name, group_map


def cubic_BCC_3d(model_name="RVE", x=0.0, y=0.0, z=0.0, dx=1.0, dy=1.0,
                  dz=1.0, lc=None, gui=False):
    # initialize gmsh and add a new model
    if not gmsh.isInitialized():
        gmsh.initialize()
    gmsh.model.add(model_name)

    # add points on corners
    point_1 = gmsh.model.geo.add_point(x - dx/2, y - dy/2, z - dz/2)
    point_2 = gmsh.model.geo.add_point(x + dx/2, y - dy/2, z - dz/2)
    point_3 = gmsh.model.geo.add_point(x + dx/2, y + dy/2, z - dz/2)
    point_4 = gmsh.model.geo.add_point(x - dx/2, y + dy/2, z - dz/2)
    point_5 = gmsh.model.geo.add_point(x - dx/2, y - dy/2, z + dz/2)
    point_6 = gmsh.model.geo.add_point(x + dx/2, y - dy/2, z + dz/2)
    point_7 = gmsh.model.geo.add_point(x + dx/2, y + dy/2, z + dz/2)
    point_8 = gmsh.model.geo.add_point(x - dx/2, y + dy/2, z + dz/2)

    # add points on surfaces
    point_9 = gmsh.model.geo.add_point(x, y, z - dz/2)
    point_10 = gmsh.model.geo.add_point(x, y, z + dz/2)
    point_11 = gmsh.model.geo.add_point(x - dx/2, y, z)
    point_12 = gmsh.model.geo.add_point(x + dx/2, y, z)

    # add lines on edges
    edge_1 = gmsh.model.geo.add_line(point_1, point_2)
    edge_2 = gmsh.model.geo.add_line(point_2, point_3)
    edge_3 = gmsh.model.geo.add_line(point_3, point_4)
    edge_4 = gmsh.model.geo.add_line(point_4, point_1)
    edge_5 = gmsh.model.geo.add_line(point_5, point_6)
    edge_6 = gmsh.model.geo.add_line(point_6, point_7)
    edge_7 = gmsh.model.geo.add_line(point_7, point_8)
    edge_8 = gmsh.model.geo.add_line(point_8, point_5)
    edge_9 = gmsh.model.geo.add_line(point_1, point_5)
    edge_10 = gmsh.model.geo.add_line(point_2, point_6)
    edge_11 = gmsh.model.geo.add_line(point_3, point_7)
    edge_12 = gmsh.model.geo.add_line(point_4, point_8)

    # add lines on surfaces
    edge_13 = gmsh.model.geo.add_line(point_1, point_9)
    edge_14 = gmsh.model.geo.add_line(point_2, point_9)
    edge_15 = gmsh.model.geo.add_line(point_3, point_9)
    edge_16 = gmsh.model.geo.add_line(point_4, point_9)
    edge_17 = gmsh.model.geo.add_line(point_5, point_10)
    edge_18 = gmsh.model.geo.add_line(point_6, point_10)
    edge_19 = gmsh.model.geo.add_line(point_7, point_10)
    edge_20 = gmsh.model.geo.add_line(point_8, point_10)
    edge_21 = gmsh.model.geo.add_line(point_5, point_11)
    edge_22 = gmsh.model.geo.add_line(point_1, point_11)
    edge_23 = gmsh.model.geo.add_line(point_4, point_11)
    edge_24 = gmsh.model.geo.add_line(point_8, point_11)
    edge_25 = gmsh.model.geo.add_line(point_6, point_12)
    edge_26 = gmsh.model.geo.add_line(point_2, point_12)
    edge_27 = gmsh.model.geo.add_line(point_3, point_12)
    edge_28 = gmsh.model.geo.add_line(point_7, point_12)

    gmsh.model.geo.synchronize()

    # create map containing information on mesh periodicity and physical groups
    group_map = {"POINT_1": (0, [point_1]),
                 "POINT_2": (0, [point_2]),
                 "POINT_3": (0, [point_3]),
                 "POINT_4": (0, [point_4]),
                 "POINT_5": (0, [point_5]),
                 "POINT_6": (0, [point_6]),
                 "POINT_7": (0, [point_7]),
                 "POINT_8": (0, [point_8]),
                 "EDGE_1": (1, [edge_1]),
                 "EDGE_2": (1, [edge_2]),
                 "EDGE_3": (1, [edge_3]),
                 "EDGE_4": (1, [edge_4]),
                 "EDGE_5": (1, [edge_5]),
                 "EDGE_6": (1, [edge_6]),
                 "EDGE_7": (1, [edge_7]),
                 "EDGE_8": (1, [edge_8]),
                 "EDGE_9": (1, [edge_9]),
                 "EDGE_10": (1, [edge_10]),
                 "EDGE_11": (1, [edge_11]),
                 "EDGE_12": (1, [edge_12]),
                 "FACE_1": (1, [edge_13, edge_14, edge_15, edge_16]),
                 "FACE_2": (1, [edge_17, edge_18, edge_19, edge_20]),
                 "FACE_3": (1, []),
                 "FACE_4": (1, []),
                 "FACE_5": (1, [edge_21, edge_22, edge_23, edge_24]),
                 "FACE_6": (1, [edge_25, edge_26, edge_27, edge_28])}

    add_boundary_groups(model_name, group_map)
    add_periodic_constraints(model_name, group_map)

    # set mesh size
    if lc is not None:
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeMin", lc)
        gmsh.option.setNumber("Mesh.MeshSizeMax", lc)

    gmsh.model.mesh.generate(1)

    if gui:
        gmsh.fltk.run()

    return model_name, group_map
