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
    group_map = {"POINT_1": [1], 
                 "POINT_2": [2], 
                 "POINT_3": [4], 
                 "POINT_4": [3],
                 "EDGE_1": [1],
                 "EDGE_2": [3],
                 "EDGE_3": [4],
                 "EDGE_4": [2]}
    
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
    gmsh.model.mesh.generate(2)
    
    if gui:
        gmsh.fltk.run()
    
    return model_name, group_map