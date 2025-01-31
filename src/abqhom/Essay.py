import gmsh
from abqhom.RVE import (add_boundary_groups, add_periodic_constraints,
                        create_boundary_sets)
from abqhom.examples import *


def model_Body_centered_cube_RVE(model_name="RVE", x=0.0, y=0.0, z=0.0, dx=1.0, dy=1.0,
                   dz=1.0, lc=None, gui=False):
    # initialize gmsh and add a new model
    if not gmsh.isInitialized():
        gmsh.initialize()
    gmsh.model.add(model_name)

    # create Geometry

    gmsh.model.occ.addBox(x-dx/2, y-dy/2, z-dz/2, dx, dy, dz,1)
    gmsh.model.occ.addLine(3,6,13)
    gmsh.model.occ.addLine(4,5,14)
    gmsh.model.occ.addLine(1,8,15)
    gmsh.model.occ.addLine(2,7,16)




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
                 "FACE_6": [6],
                 }

    add_boundary_groups(model_name, group_map)
    add_periodic_constraints(model_name, group_map, dx, dy, dz)

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

def model_face_centered_cube_RVE(model_name="RVE", x=0.0, y=0.0, z=0.0, dx=1.0, dy=1.0,
                   dz=1.0, lc=None, gui=False):
    if not gmsh.isInitialized():
        gmsh.initialize()
    gmsh.model.add(model_name)

    gmsh.model.occ.addBox(x-dx/2, y-dy/2, z-dz/2, dx, dy, dz,1)
    gmsh.model.occ.addPoint(x+dx/2,y,z,lc,9)
    gmsh.model.occ.addPoint(x-dx/2,y,z,lc,10)
    gmsh.model.occ.addPoint(x,y+dy/2,z,lc,11)
    gmsh.model.occ.addPoint(x,y-dy/2,z,lc,12)
    gmsh.model.occ.addPoint(x,y,z+dz/2,lc,13)
    gmsh.model.occ.addPoint(x,y,z-dz/2,lc,14)
    gmsh.model.occ.fuse([(0, 9)], [(2,2)])
    gmsh.model.occ.fuse([(0, 10)], [(2,1)])
    gmsh.model.occ.fuse([(0, 11)], [(2,4)])
    gmsh.model.occ.fuse([(0, 12)], [(2,3)])
    gmsh.model.occ.fuse([(0, 13)], [(2,6)])
    gmsh.model.occ.fuse([(0, 14)], [(2,5)])
    gmsh.model.occ.addLine(9,10,13)
    gmsh.model.occ.addLine(11,12,14)
    gmsh.model.occ.addLine(13,14,15)
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
                 "FACE_6": [6],
                 }
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

def model_all_face_centered_cubic(model_name="RVE", x=0.0, y=0.0, z=0.0, dx=1.0, dy=1.0,
                   dz=1.0, lc=None, gui=False):
    if not gmsh.isInitialized():
        gmsh.initialize()
    gmsh.model.add(model_name)

    gmsh.model.occ.addBox(x-dx/2, y-dy/2, z-dz/2, dx, dy, dz,1)
    gmsh.model.occ.addLine(1,6,13)
    gmsh.model.occ.addLine(1,4,14)
    gmsh.model.occ.addLine(1,7,15)
    gmsh.model.occ.addLine(2,3,16)
    gmsh.model.occ.addLine(2,5,17)
    gmsh.model.occ.addLine(2,8,18)
    gmsh.model.occ.addLine(3,5,19)
    gmsh.model.occ.addLine(3,8,20)
    gmsh.model.occ.addLine(4,6,21)
    gmsh.model.occ.addLine(4,7,22)
    gmsh.model.occ.addLine(5,8,23)
    gmsh.model.occ.addLine(6,7,24)
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
                 "FACE_6": [6],
                 }
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

def model_diamond_lattice_structure(model_name="RVE", x=0.0, y=0.0, z=0.0, dx=1.0, dy=1.0,
                   dz=1.0, lc=None, gui=False):
    if not gmsh.isInitialized():
        gmsh.initialize()
    gmsh.model.add(model_name)

    gmsh.model.occ.addBox(x - dx / 2, y - dy / 2, z - dz / 2, dx, dy, dz, 1)
    gmsh.model.occ.addPoint(x + dx / 2, y, z, lc, 9)
    gmsh.model.occ.addPoint(x - dx / 2, y, z, lc, 10)
    gmsh.model.occ.addPoint(x, y + dy / 2, z, lc, 11)
    gmsh.model.occ.addPoint(x, y - dy / 2, z, lc, 12)
    gmsh.model.occ.addPoint(x, y, z + dz / 2, lc, 13)
    gmsh.model.occ.addPoint(x, y, z - dz / 2, lc, 14)
    gmsh.model.occ.addLine(9,11,13)
    gmsh.model.occ.addLine(9,12,14)
    gmsh.model.occ.addLine(10,11,15)
    gmsh.model.occ.addLine(10,12,16)
    gmsh.model.occ.addLine(11,13,17)
    gmsh.model.occ.addLine(11,14,18)
    gmsh.model.occ.addLine(12,13,19)
    gmsh.model.occ.addLine(12,14,20)
    gmsh.model.occ.addLine(11,12,21)
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
                 "FACE_6": [6],
                 }
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
def model_ACC_lattice_structure(model_name="RVE", x=0.0, y=0.0, z=0.0, dx=1.0, dy=1.0,
                   dz=1.0, lc=None, gui=False):
    if not gmsh.isInitialized():
        gmsh.initialize()
    gmsh.model.add(model_name)

    gmsh.model.occ.addBox(x - dx / 2, y - dy / 2, z - dz / 2, dx, dy, dz, 1)
    gmsh.model.occ.addPoint(x , y+dy/2, z+dz/2, lc, 9)
    gmsh.model.occ.addPoint(x , y-dy/2, z-dz/2, lc, 10)
    gmsh.model.occ.addPoint(x, y+dy/2, z-dz/2, lc, 11)
    gmsh.model.occ.addPoint(x, y-dy/2, z+dz/2, lc, 12)
    gmsh.model.occ.addPoint(x+dx/2, y, z+dz/2, lc, 13)
    gmsh.model.occ.addPoint(x-dx/2, y, z-dz/2, lc, 14)
    gmsh.model.occ.addPoint(x+dx/2, y, z-dz/2, lc, 15)
    gmsh.model.occ.addPoint(x-dx/2, y, z+dz/2, lc, 16)
    gmsh.model.occ.addPoint(x+dx/2, y+dy/2, z, lc, 17)
    gmsh.model.occ.addPoint(x-dx/2, y-dy/2, z, lc, 18)
    gmsh.model.occ.addPoint(x+dx/2, y-dy/2, z, lc, 19)
    gmsh.model.occ.addPoint(x-dx/2, y+dy/2, z, lc, 20)
    gmsh.model.occ.addLine(9,10,13)
    gmsh.model.occ.addLine(11,12,14)
    gmsh.model.occ.addLine(13,14,15)
    gmsh.model.occ.addLine(15,16,16)
    gmsh.model.occ.addLine(17,18,17)
    gmsh.model.occ.addLine(19,20,18)


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
                 "FACE_6": [6],
                 }
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

def model_central_cube_lattice_structure(model_name="RVE", x=0.0, y=0.0, z=0.0, dx=1.0, dy=1.0,
                   dz=1.0,ratio_x=1/2,ratio_y=1/2,ratio_z=1/2, lc=None, gui=False):
    if not gmsh.isInitialized():
        gmsh.initialize()
    gmsh.model.add(model_name)

    gmsh.model.occ.addBox(x - dx / 2, y - dy / 2, z - dz / 2, dx, dy, dz, 1)
    dx1,dy1,dz1=dx*ratio_x,dy*ratio_y,dz*ratio_z
    gmsh.model.occ.addBox(x - dx1 / 2, y - dy1 / 2, z - dz1 / 2, dx1/2, dy1/2, dz1/2, 2)
    gmsh.model.occ.addLine(1,9)
    gmsh.model.occ.addLine(2,10)
    gmsh.model.occ.addLine(3,11)
    gmsh.model.occ.addLine(4,12)
    gmsh.model.occ.addLine(5,13)
    gmsh.model.occ.addLine(6,14)
    gmsh.model.occ.addLine(7,15)
    gmsh.model.occ.addLine(8,16)

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
                 "FACE_6": [6],
                 }
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


dx = dy = dz = 10
lc = 5
model_Body_centered_cube_RVE("RVE", 0.0, 0.0, 0.0, dx, dy,dz,lc, True)
