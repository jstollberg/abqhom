import gmsh
from abqhom.RVE import (add_boundary_groups, add_periodic_constraints,
                        create_boundary_sets)
from abqhom.examples import *


def model_Body_centered_cube_RVE(model_name="BCC", x=0.0, y=0.0, z=0.0, dx=1.0, dy=1.0,
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

    gmsh.write("{}.geo".format(model_name))

    if gui:
        gmsh.fltk.run()

    return model_name, group_map

def model_face_centered_cube_RVE(model_name="FCC", x=0.0, y=0.0, z=0.0, dx=1.0, dy=1.0,
                   dz=1.0, lc=None, gui=False):
    if not gmsh.isInitialized():
        gmsh.initialize()
    gmsh.model.add(model_name)


    Dx, Dy, Dz = dx / 2, dy / 2, dz / 2

    gmsh.model.occ.addPoint(x - Dx, y - Dy, z + Dz, lc, 1)
    gmsh.model.occ.addPoint(x - Dx, y - Dy, z - Dz, lc, 2)
    gmsh.model.occ.addPoint(x - Dx, y + Dy, z - Dz, lc, 3)
    gmsh.model.occ.addPoint(x - Dx, y + Dy, z + Dz, lc, 4)
    gmsh.model.occ.addPoint(x + Dx, y - Dy, z + Dz, lc, 5)
    gmsh.model.occ.addPoint(x + Dx, y - Dy, z - Dz, lc, 6)
    gmsh.model.occ.addPoint(x + Dx, y + Dy, z - Dz, lc, 7)
    gmsh.model.occ.addPoint(x + Dx, y + Dy, z + Dz, lc, 8)

    gmsh.model.occ.addPoint(x+Dx,y,z,lc,9)
    gmsh.model.occ.addPoint(x-Dx,y,z,lc,10)
    gmsh.model.occ.addPoint(x,y+Dy,z,lc,11)
    gmsh.model.occ.addPoint(x,y-Dy,z,lc,12)
    gmsh.model.occ.addPoint(x,y,z+Dz,lc,13)
    gmsh.model.occ.addPoint(x,y,z-Dz,lc,14)

    gmsh.model.occ.addLine(1,2, 1)
    gmsh.model.occ.addLine(2,3, 2)
    gmsh.model.occ.addLine(3,4, 3)
    gmsh.model.occ.addLine(4,1, 4)
    gmsh.model.occ.addLine(5,6, 5)
    gmsh.model.occ.addLine(6, 7, 6)
    gmsh.model.occ.addLine(7,  8, 7)
    gmsh.model.occ.addLine(8, 5, 8)
    gmsh.model.occ.addLine(2, 6, 9)
    gmsh.model.occ.addLine(1, 5, 10)
    gmsh.model.occ.addLine(3, 7, 11)
    gmsh.model.occ.addLine(4,8, 12)

    gmsh.model.occ.addLine(9,10,13)
    gmsh.model.occ.addLine(11,12,14)
    gmsh.model.occ.addLine(13,14,15)

    gmsh.model.occ.addCurveLoop([1, 2, 3, 4], 1)
    gmsh.model.occ.addCurveLoop([5, 6, 7, 8], 2)
    gmsh.model.occ.addCurveLoop([10, 5, 9, 1], 3)
    gmsh.model.occ.addCurveLoop([12, 7, 11, 3], 4)
    gmsh.model.occ.addCurveLoop([9, 6, 11, 2], 5)
    gmsh.model.occ.addCurveLoop([10, 8, 12, 4], 6)



    gmsh.model.occ.addSurfaceFilling(1, 1,[10])
    gmsh.model.occ.addSurfaceFilling(2, 2,[9])
    gmsh.model.occ.addSurfaceFilling(3, 3,[12])
    gmsh.model.occ.addSurfaceFilling(4, 4,[11])
    gmsh.model.occ.addSurfaceFilling(5, 5,[14])
    gmsh.model.occ.addSurfaceFilling(6, 6,[13])

    gmsh.model.occ.addSurfaceLoop([1, 2, 3, 4, 5, 6], 1)
    gmsh.model.occ.addVolume([1], 1)

    gmsh.model.addPhysicalGroup(0, [9])
    gmsh.model.addPhysicalGroup(0, [10])
    gmsh.model.addPhysicalGroup(0, [11])
    gmsh.model.addPhysicalGroup(0, [12])
    gmsh.model.addPhysicalGroup(0, [13])
    gmsh.model.addPhysicalGroup(0, [14])

    gmsh.model.addPhysicalGroup(1, [13])
    gmsh.model.addPhysicalGroup(1, [14])
    gmsh.model.addPhysicalGroup(1, [15])





    gmsh.model.occ.synchronize()

    # create map containing information on mesh periodicity and physical groups
    group_map = {"POINT_1": [1],
                 "POINT_2": [2],
                 "POINT_3": [3],
                 "POINT_4": [4],
                 "POINT_5": [5],
                 "POINT_6": [6],
                 "POINT_7": [7],
                 "POINT_8": [8],
                 "EDGE_1": [1],
                 "EDGE_2": [2],
                 "EDGE_3": [3],
                 "EDGE_4": [4],
                 "EDGE_5": [5],
                 "EDGE_6": [6],
                 "EDGE_7": [7],
                 "EDGE_8": [8],
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

    gmsh.write("{}.geo".format(model_name))

    if gui:
        gmsh.fltk.run()

    return model_name, group_map

def model_all_face_centered_cubic(model_name="AFCC", x=0.0, y=0.0, z=0.0, dx=1.0, dy=1.0,
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

    gmsh.write("{}.geo_unrolled".format(model_name))

    if gui:
        gmsh.fltk.run()

    return model_name, group_map

def model_diamond_lattice_structure(model_name="Diamond", x=0.0, y=0.0, z=0.0, dx=1.0, dy=1.0,
                   dz=1.0, lc=None, gui=False):
    if not gmsh.isInitialized():
        gmsh.initialize()
    gmsh.model.add(model_name)

    Dx, Dy, Dz = dx / 2, dy / 2, dz / 2

    gmsh.model.occ.addPoint(x - Dx, y - Dy, z + Dz, lc, 1)
    gmsh.model.occ.addPoint(x - Dx, y - Dy, z - Dz, lc, 2)
    gmsh.model.occ.addPoint(x - Dx, y + Dy, z - Dz, lc, 3)
    gmsh.model.occ.addPoint(x - Dx, y + Dy, z + Dz, lc, 4)
    gmsh.model.occ.addPoint(x + Dx, y - Dy, z + Dz, lc, 5)
    gmsh.model.occ.addPoint(x + Dx, y - Dy, z - Dz, lc, 6)
    gmsh.model.occ.addPoint(x + Dx, y + Dy, z - Dz, lc, 7)
    gmsh.model.occ.addPoint(x + Dx, y + Dy, z + Dz, lc, 8)

    gmsh.model.occ.addLine(1, 2, 1)
    gmsh.model.occ.addLine(2, 3, 2)
    gmsh.model.occ.addLine(3, 4, 3)
    gmsh.model.occ.addLine(4, 1, 4)
    gmsh.model.occ.addLine(5, 6, 5)
    gmsh.model.occ.addLine(6, 7, 6)
    gmsh.model.occ.addLine(7, 8, 7)
    gmsh.model.occ.addLine(8, 5, 8)
    gmsh.model.occ.addLine(2, 6, 9)
    gmsh.model.occ.addLine(1, 5, 10)
    gmsh.model.occ.addLine(3, 7, 11)
    gmsh.model.occ.addLine(4, 8, 12)

    gmsh.model.occ.addPoint(x + dx / 2, y, z,lc, 9)
    gmsh.model.occ.addPoint(x - dx / 2, y, z,lc, 10)
    gmsh.model.occ.addPoint(x, y + dy / 2, z,lc, 11)
    gmsh.model.occ.addPoint(x, y - dy / 2, z,lc, 12)
    gmsh.model.occ.addPoint(x, y, z + dz / 2,lc, 13)
    gmsh.model.occ.addPoint(x, y, z - dz / 2,lc, 14)


    gmsh.model.occ.addLine(9,11,13)
    gmsh.model.occ.addLine(9,12,14)
    gmsh.model.occ.addLine(10,11,15)
    gmsh.model.occ.addLine(10,12,16)
    gmsh.model.occ.addLine(11,13,17)
    gmsh.model.occ.addLine(11,14,18)
    gmsh.model.occ.addLine(12,13,19)
    gmsh.model.occ.addLine(12,14,20)
    gmsh.model.occ.addLine(11,12,21)



    gmsh.model.occ.addCurveLoop([1, 2, 3, 4],1)
    gmsh.model.occ.addCurveLoop([5, 6, 7, 8], 2)
    gmsh.model.occ.addCurveLoop([10, 5, 9, 1], 3)
    gmsh.model.occ.addCurveLoop([12, 7, 11, 3], 4)
    gmsh.model.occ.addCurveLoop([9, 6, 11, 2], 5)
    gmsh.model.occ.addCurveLoop([10, 8, 12, 4], 6)

    gmsh.model.occ.addSurfaceFilling(1, 1, [10])
    gmsh.model.occ.addSurfaceFilling(2, 2, [9])
    gmsh.model.occ.addSurfaceFilling(3, 3, [12])
    gmsh.model.occ.addSurfaceFilling(4, 4, [11])
    gmsh.model.occ.addSurfaceFilling(5, 5, [14])
    gmsh.model.occ.addSurfaceFilling(6, 6, [13])


    gmsh.model.addPhysicalGroup(0, [9])
    gmsh.model.addPhysicalGroup(0, [10])
    gmsh.model.addPhysicalGroup(0, [11])
    gmsh.model.addPhysicalGroup(0, [12])
    gmsh.model.addPhysicalGroup(0, [13])
    gmsh.model.addPhysicalGroup(0, [14])

    gmsh.model.addPhysicalGroup(1, [13])
    gmsh.model.addPhysicalGroup(1, [14])
    gmsh.model.addPhysicalGroup(1, [15])
    gmsh.model.addPhysicalGroup(1, [16])
    gmsh.model.addPhysicalGroup(1, [17])
    gmsh.model.addPhysicalGroup(1, [18])
    gmsh.model.addPhysicalGroup(1, [19])
    gmsh.model.addPhysicalGroup(1, [20])
    gmsh.model.addPhysicalGroup(1, [21])

    gmsh.model.occ.addSurfaceLoop([1, 2, 3, 4, 5, 6], 1)
    gmsh.model.occ.addVolume([1], 1)




    gmsh.model.occ.synchronize()


    # create map containing information on mesh periodicity and physical groups
    group_map = {"POINT_1": [1],
                 "POINT_2": [2],
                 "POINT_3": [3],
                 "POINT_4": [4],
                 "POINT_5": [5],
                 "POINT_6": [6],
                 "POINT_7": [7],
                 "POINT_8": [8],
                 "EDGE_1": [1],
                 "EDGE_2": [2],
                 "EDGE_3": [3],
                 "EDGE_4": [4],
                 "EDGE_5": [5],
                 "EDGE_6": [6],
                 "EDGE_7": [7],
                 "EDGE_8": [8],
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

    gmsh.write("{}.msh".format(model_name))

    if gui:
        gmsh.fltk.run()

    return model_name, group_map
def model_ACC_lattice_structure(model_name="ACC", x=0.0, y=0.0, z=0.0, dx=1.0, dy=1.0,
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
def model_ACC_lattice_structure1(model_name="ACC", x=0.0, y=0.0, z=0.0, dx=1.0, dy=1.0,
                   dz=1.0, lc=None, gui=False):
    if not gmsh.isInitialized():
        gmsh.initialize()
    gmsh.model.add(model_name)
    Dx,Dy,Dz=dx/2,dy/2,dz/2


    gmsh.model.occ.addPoint(x-Dx,y-Dy,z+Dz,lc, 1)
    gmsh.model.occ.addPoint(x-Dx,y-Dy,z-Dz,lc, 2)
    gmsh.model.occ.addPoint(x-Dx,y+Dy,z-Dz,lc, 3)
    gmsh.model.occ.addPoint(x-Dx,y+Dy,z+Dz,lc, 4)
    gmsh.model.occ.addPoint(x+Dx,y-Dy,z+Dz,lc, 5)
    gmsh.model.occ.addPoint(x+Dx,y-Dy,z-Dz,lc, 6)
    gmsh.model.occ.addPoint(x+Dx,y+Dy,z-Dz,lc, 7)

    gmsh.model.occ.addPoint(x+Dx,y+Dy,z+Dz,lc, 8)
    gmsh.model.occ.addPoint(x-Dx,y-Dy,z,lc, 9)
    gmsh.model.occ.addPoint(x-Dx,y,z-Dz,lc, 10)
    gmsh.model.occ.addPoint(x-Dx,y+Dy,z,lc, 11)
    gmsh.model.occ.addPoint(x-Dx,y,z+Dz,lc, 12)
    gmsh.model.occ.addPoint(x+Dx,y-Dy,z,lc, 13)
    gmsh.model.occ.addPoint(x+Dx,y,z-Dz,lc, 14)
    gmsh.model.occ.addPoint(x+Dx,y+Dy,z,lc, 15)
    gmsh.model.occ.addPoint(x+Dx,y,z+Dz,lc, 16)
    gmsh.model.occ.addPoint(x,y-Dy,z-Dz,lc, 17)
    gmsh.model.occ.addPoint(x,y-Dy,z+Dz,lc, 18)
    gmsh.model.occ.addPoint(x,y+Dy,z-Dz,lc, 19)
    gmsh.model.occ.addPoint(x,y+Dy,z+Dz,lc, 20)

    gmsh.model.occ.addLine(1, 9, 1)
    gmsh.model.occ.addLine(2, 10, 2)
    gmsh.model.occ.addLine(3, 11, 3)
    gmsh.model.occ.addLine(4, 12, 4)
    gmsh.model.occ.addLine(5, 13, 5)
    gmsh.model.occ.addLine(6, 14, 6)
    gmsh.model.occ.addLine(7, 15, 7)
    gmsh.model.occ.addLine(8, 16, 8)
    gmsh.model.occ.addLine(2, 17, 9)
    gmsh.model.occ.addLine(1, 18, 10)
    gmsh.model.occ.addLine(3, 19, 11)
    gmsh.model.occ.addLine(4, 20, 12)
    gmsh.model.occ.addLine(9, 2, 13)
    gmsh.model.occ.addLine(10, 3, 14)
    gmsh.model.occ.addLine(11, 4, 15)
    gmsh.model.occ.addLine(12, 1, 16)
    gmsh.model.occ.addLine(13, 6, 17)
    gmsh.model.occ.addLine(14, 7, 18)
    gmsh.model.occ.addLine(15, 8, 19)
    gmsh.model.occ.addLine(16, 5, 20)
    gmsh.model.occ.addLine(17, 6, 21)
    gmsh.model.occ.addLine(18, 5, 22)
    gmsh.model.occ.addLine(19, 7, 23)
    gmsh.model.occ.addLine(20, 8, 24)

    gmsh.model.occ.addLine(9, 15, 25)
    gmsh.model.occ.addLine(10, 16, 26)
    gmsh.model.occ.addLine(11, 13, 27)
    gmsh.model.occ.addLine(12, 14, 28)
    gmsh.model.occ.addLine(17, 20, 29)
    gmsh.model.occ.addLine(18, 19, 30)




    gmsh.model.occ.addCurveLoop([1,13, 2,14, 3,15, 4,16], 1)
    gmsh.model.occ.addCurveLoop([5,17, 6,18, 7,19, 8,20], 2)
    gmsh.model.occ.addCurveLoop([10,22, 5,17, 21,9, 13,1], 3)
    gmsh.model.occ.addCurveLoop([12,24, 19,7, 23,11, 3,15], 4)
    gmsh.model.occ.addCurveLoop([9,21, 6,18, 23,11, 14,2], 5)
    gmsh.model.occ.addCurveLoop([10,22, 20,8, 24,12, 4,16], 6)

    gmsh.model.occ.addPlaneSurface([1], 1)
    gmsh.model.occ.addPlaneSurface([2], 2)
    gmsh.model.occ.addPlaneSurface([3], 3)
    gmsh.model.occ.addPlaneSurface([4], 4)
    gmsh.model.occ.addPlaneSurface([5], 5)
    gmsh.model.occ.addPlaneSurface([6], 6)





    gmsh.model.occ.addSurfaceLoop([1,2,3,4,5,6],1)
    gmsh.model.occ.addVolume([1],1)

    gmsh.model.addPhysicalGroup(1, [25])
    gmsh.model.addPhysicalGroup(1, [26])
    gmsh.model.addPhysicalGroup(1, [27])
    gmsh.model.addPhysicalGroup(1, [28])
    gmsh.model.addPhysicalGroup(1, [29])
    gmsh.model.addPhysicalGroup(1, [30])

    gmsh.model.addPhysicalGroup(0, [9])
    gmsh.model.addPhysicalGroup(0, [10])
    gmsh.model.addPhysicalGroup(0, [11])
    gmsh.model.addPhysicalGroup(0, [12])
    gmsh.model.addPhysicalGroup(0, [13])
    gmsh.model.addPhysicalGroup(0, [14])
    gmsh.model.addPhysicalGroup(0, [15])
    gmsh.model.addPhysicalGroup(0, [16])
    gmsh.model.addPhysicalGroup(0, [17])
    gmsh.model.addPhysicalGroup(0, [18])
    gmsh.model.addPhysicalGroup(0, [19])
    gmsh.model.addPhysicalGroup(0, [20])

    gmsh.model.occ.synchronize()

    # create map containing information on mesh periodicity and physical groups
    group_map = {"POINT_1": [1],
                 "POINT_2": [2],
                 "POINT_3": [3],
                 "POINT_4": [4],
                 "POINT_5": [5],
                 "POINT_6": [6],
                 "POINT_7": [7],
                 "POINT_8": [8],
                 "EDGE_1": [1,13],
                 "EDGE_2": [2,14],
                 "EDGE_3": [3,15],
                 "EDGE_4": [4,16],
                 "EDGE_5": [5,17],
                 "EDGE_6": [6,18],
                 "EDGE_7": [7,19],
                 "EDGE_8": [8,20],
                 "EDGE_9": [9,21],
                 "EDGE_10": [10,22],
                 "EDGE_11": [11,23],
                 "EDGE_12": [12,24],
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

    gmsh.write("{}.msh".format(model_name))

    if gui:
        gmsh.fltk.run()

    return model_name, group_map


def model_ACC_lattice_structure2(model_name="ACC", x=0.0, y=0.0, z=0.0, dx=1.0, dy=1.0,
                   dz=1.0, lc=None, gui=False):
    if not gmsh.isInitialized():
        gmsh.initialize()
    gmsh.model.add(model_name)
    Dx,Dy,Dz=dx/2,dy/2,dz/2


    gmsh.model.occ.addPoint(x-Dx,y-Dy,z+Dz,lc, 1)
    gmsh.model.occ.addPoint(x-Dx,y-Dy,z-Dz,lc, 2)
    gmsh.model.occ.addPoint(x-Dx,y+Dy,z-Dz,lc, 3)
    gmsh.model.occ.addPoint(x-Dx,y+Dy,z+Dz,lc, 4)
    gmsh.model.occ.addPoint(x+Dx,y-Dy,z+Dz,lc, 5)
    gmsh.model.occ.addPoint(x+Dx,y-Dy,z-Dz,lc, 6)
    gmsh.model.occ.addPoint(x+Dx,y+Dy,z-Dz,lc, 7)
    gmsh.model.occ.addPoint(x+Dx,y+Dy,z+Dz,lc, 8)



    gmsh.model.occ.addPoint(x-Dx,y-Dy,z,lc, 9)
    gmsh.model.occ.addPoint(x-Dx,y,z-Dz,lc, 10)
    gmsh.model.occ.addPoint(x-Dx,y+Dy,z,lc, 11)
    gmsh.model.occ.addPoint(x-Dx,y,z+Dz,lc, 12)
    gmsh.model.occ.addPoint(x+Dx,y-Dy,z,lc, 13)
    gmsh.model.occ.addPoint(x+Dx,y,z-Dz,lc, 14)
    gmsh.model.occ.addPoint(x+Dx,y+Dy,z,lc, 15)
    gmsh.model.occ.addPoint(x+Dx,y,z+Dz,lc, 16)
    gmsh.model.occ.addPoint(x,y-Dy,z-Dz,lc, 17)
    gmsh.model.occ.addPoint(x,y-Dy,z+Dz,lc, 18)
    gmsh.model.occ.addPoint(x,y+Dy,z-Dz,lc, 19)
    gmsh.model.occ.addPoint(x,y+Dy,z+Dz,lc, 20)

    gmsh.model.occ.addSpline([1, 9, 2], 1)
    gmsh.model.occ.addSpline([2, 10, 3], 2)
    gmsh.model.occ.addSpline([3, 11, 4], 3)
    gmsh.model.occ.addSpline([4, 12, 1], 4)
    gmsh.model.occ.addSpline([5, 13, 6], 5)
    gmsh.model.occ.addSpline([6, 14, 7], 6)
    gmsh.model.occ.addSpline([7, 15, 8], 7)
    gmsh.model.occ.addSpline([8, 16, 5], 8)
    gmsh.model.occ.addSpline([2, 17, 6], 9)
    gmsh.model.occ.addSpline([1, 18, 5], 10)
    gmsh.model.occ.addSpline([3, 19, 7], 11)
    gmsh.model.occ.addSpline([4, 20, 8], 12)

    gmsh.model.occ.addLine(9, 15, 13)
    gmsh.model.occ.addLine(10, 16, 14)
    gmsh.model.occ.addLine(11, 13, 15)
    gmsh.model.occ.addLine(12, 14, 16)
    gmsh.model.occ.addLine(17, 20, 17)
    gmsh.model.occ.addLine(18, 19, 18)



    gmsh.model.occ.addCurveLoop([1, 2, 3, 4], 1)
    gmsh.model.occ.addCurveLoop([5, 6, 7, 8], 2)
    gmsh.model.occ.addCurveLoop([10, 5, 9, 1], 3)
    gmsh.model.occ.addCurveLoop([12, 7, 11, 3], 4)
    gmsh.model.occ.addCurveLoop([9, 6, 11, 2], 5)
    gmsh.model.occ.addCurveLoop([10, 8, 12, 4], 6)



    gmsh.model.occ.addSurfaceFilling(1, 1)
    gmsh.model.occ.addSurfaceFilling(2, 2)
    gmsh.model.occ.addSurfaceFilling(3, 3)
    gmsh.model.occ.addSurfaceFilling(4, 4)
    gmsh.model.occ.addSurfaceFilling(5, 5)
    gmsh.model.occ.addSurfaceFilling(6, 6)



    gmsh.model.occ.addSurfaceLoop([1,2,3,4,5,6],1)
    gmsh.model.occ.addVolume([1],1)

    gmsh.model.occ.fragment([(1, 1)], [(1, 13)])
    gmsh.model.occ.fragment([(1, 2)], [(1, 14)])
    gmsh.model.occ.fragment([(1, 3)], [(1, 15)])
    gmsh.model.occ.fragment([(1, 4)], [(1, 16)])
    gmsh.model.occ.fragment([(1, 5)], [(1, 15)])
    gmsh.model.occ.fragment([(1, 6)], [(1, 16)])
    gmsh.model.occ.fragment([(1, 7)], [(1, 13)])
    gmsh.model.occ.fragment([(1, 8)], [(1, 14)])
    gmsh.model.occ.fragment([(1, 9)], [(1, 17)])
    gmsh.model.occ.fragment([(1, 10)], [(1, 18)])
    gmsh.model.occ.fragment([(1, 11)], [(1, 18)])
    gmsh.model.occ.fragment([(1, 12)], [(1, 17)])


    gmsh.model.addPhysicalGroup(1,[13])
    gmsh.model.addPhysicalGroup(1,[14])
    gmsh.model.addPhysicalGroup(1,[15])
    gmsh.model.addPhysicalGroup(1,[16])
    gmsh.model.addPhysicalGroup(1,[17])
    gmsh.model.addPhysicalGroup(1,[18])


    gmsh.model.addPhysicalGroup(0,[9])
    gmsh.model.addPhysicalGroup(0, [10])
    gmsh.model.addPhysicalGroup(0, [11])
    gmsh.model.addPhysicalGroup(0, [12])
    gmsh.model.addPhysicalGroup(0, [13])
    gmsh.model.addPhysicalGroup(0, [14])
    gmsh.model.addPhysicalGroup(0,[15])
    gmsh.model.addPhysicalGroup(0,[16])
    gmsh.model.addPhysicalGroup(0,[17])
    gmsh.model.addPhysicalGroup(0,[18])
    gmsh.model.addPhysicalGroup(0,[19])
    gmsh.model.addPhysicalGroup(0,[20])



    gmsh.model.occ.synchronize()

    # create map containing information on mesh periodicity and physical groups
    group_map = {"POINT_1": [1],
                 "POINT_2": [2],
                 "POINT_3": [3],
                 "POINT_4": [4],
                 "POINT_5": [5],
                 "POINT_6": [6],
                 "POINT_7": [7],
                 "POINT_8": [8],
                 "EDGE_1": [1],
                 "EDGE_2": [2],
                 "EDGE_3": [3],
                 "EDGE_4": [4],
                 "EDGE_5": [5],
                 "EDGE_6": [6],
                 "EDGE_7": [7],
                 "EDGE_8": [8],
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

def model_central_cube_lattice_structure(model_name="CCube", x=0.0, y=0.0, z=0.0, dx=1.0, dy=1.0,
                   dz=1.0,ratio_x=1/2,ratio_y=1/2,ratio_z=1/2, lc=None, gui=False):
    if not gmsh.isInitialized():
        gmsh.initialize()
    gmsh.model.add(model_name)

    gmsh.model.occ.addBox(x - dx / 2, y - dy / 2, z - dz / 2, dx, dy, dz, 1)
    dxa,dya,dza= dx*ratio_x, dy*ratio_y, dz*ratio_z
    gmsh.model.occ.addBox(x-dxa/2, y-dya/2, z-dza/2, dxa, dya, dza, 2)
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

    gmsh.write("{}.msh".format(model_name))

    if gui:
        gmsh.fltk.run()

    return model_name, group_map


def model_essay(model_name="RVE", x=0.0, y=0.0, z=0.0, dx=1.0, dy=1.0,
                   dz=1.0, lc=None, gui=False):
    if not gmsh.isInitialized():
        gmsh.initialize()
    gmsh.model.add(model_name)
    gmsh.model.occ.addPoint(x,y,z,lc,1)
    gmsh.model.occ.addPoint(x+dx,y,z,lc,2)
    gmsh.model.occ.addPoint(x+dx,y+dy,z,lc,3)
    gmsh.model.occ.addPoint(x,y+dy,z,lc,4)
    gmsh.model.occ.addPoint(x+dx/2,y+dy/2,z,lc,5)
    gmsh.model.occ.addPoint(x+dx/2,y-dy,z,lc,6)
    gmsh.model.occ.addSpline([1,6,2],1)
    gmsh.model.occ.addLine(2,3,2)
    gmsh.model.occ.addLine(3,4,3)
    gmsh.model.occ.addLine(4,1,4)
    gmsh.model.occ.addCurveLoop([1,2,3,4],1)
    gmsh.model.occ.addSurfaceFilling(1,1,[5])
    gmsh.model.occ.synchronize()

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

    return;


def model_CC_latttice_structure(model_name="CC", x=0.0, y=0.0, z=0.0, dx=1.0, dy=1.0,
                   dz=1.0,radius=0.25, lc=None, gui=False):
    if not gmsh.isInitialized():
        gmsh.initialize()
    gmsh.model.add(model_name)
    gmsh.model.occ.addCylinder(x-dx/2,y,z,dx,0,0,radius,1)
    gmsh.model.occ.addCylinder(x,y-dy/2,z,0,dy,0,radius,2)
    gmsh.model.occ.addCylinder(x,y,z-dz/2,0,0,dz,radius,3)
    gmsh.model.occ.fuse([(3,1)],[(3,2)],4)
    gmsh.model.occ.fuse([(3,4)],[(3,3)],5)
    gmsh.model.occ.synchronize()
    group_map = {"POINT_1": [22],
                 "POINT_2": [9],
                 "POINT_3": [21],
                 "POINT_4": [18],
                 "POINT_5": [20],
                 "POINT_6": [10],
                 "EDGE_1": [32],
                 "EDGE_2": [10],
                 "EDGE_3": [30],
                 "EDGE_4": [23],
                 "EDGE_5": [27],
                 "EDGE_6": [18],
                 "FACE_1": [12],
                 "FACE_2": [6],
                 "FACE_3": [11],
                 "FACE_4": [9],
                 "FACE_5": [10],
                 "FACE_6": [8],
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

    gmsh.write("{}.msh".format(model_name))

    if gui:
        gmsh.fltk.run()
    return model_name, group_map

def model_Octahedron_lattice_structure(model_name="Octahedron", x=0.0, y=0.0, z=0.0, dx=1.0, dy=1.0,
                   dz=1.0, lc=None, gui=False):
    if not gmsh.isInitialized():
        gmsh.initialize()
    gmsh.model.add(model_name)

    Dx,Dy,Dz=dx/2,dy/2,dz/2

    gmsh.model.occ.addPoint(x-Dx,y,z,lc,1)
    gmsh.model.occ.addPoint(x+Dx,y,z,lc,2)
    gmsh.model.occ.addPoint(x,y-Dy,z,lc,3)
    gmsh.model.occ.addPoint(x,y+Dy,z,lc,4)
    gmsh.model.occ.addPoint(x,y,z-Dz,lc,5)
    gmsh.model.occ.addPoint(x,y,z+Dz,lc,6)

    gmsh.model.occ.addLine(1,3,1)
    gmsh.model.occ.addLine(1,4,2)
    gmsh.model.occ.addLine(1,5,3)
    gmsh.model.occ.addLine(1,6,4)
    gmsh.model.occ.addLine(2,3,5)
    gmsh.model.occ.addLine(2,4,6)
    gmsh.model.occ.addLine(2,5,7)
    gmsh.model.occ.addLine(2,6,8)
    gmsh.model.occ.addLine(3,5,9)
    gmsh.model.occ.addLine(3,6,10)
    gmsh.model.occ.addLine(4,5,11)
    gmsh.model.occ.addLine(4,6,12)


    gmsh.model.occ.synchronize()

    # set mesh size
    if lc is not None:
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeMin", lc)
        gmsh.option.setNumber("Mesh.MeshSizeMax", lc)

    # create mesh using Delauny algorithm
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)
    gmsh.model.mesh.generate(3)

    gmsh.write("{}.msh".format(model_name))

    if gui:
        gmsh.fltk.run()



dx = dy = dz = 10
lc = 5
model_central_cube_lattice_structure(x=0.0, y=0.0, z=0.0, dx=dx, dy=dy,dz=dz,lc=lc, gui=True)

"""
gmsh.model.occ.addLine(1, 9, 1)
gmsh.model.occ.addLine(2, 10, 2)
gmsh.model.occ.addLine(3, 11, 3)
gmsh.model.occ.addLine(4, 12, 4)
gmsh.model.occ.addLine(5, 13, 5)
gmsh.model.occ.addLine(6, 14, 6)
gmsh.model.occ.addLine(7, 15, 7)
gmsh.model.occ.addLine(8, 16, 8)
gmsh.model.occ.addLine(2, 17, 9)
gmsh.model.occ.addLine(1, 18, 10)
gmsh.model.occ.addLine(3, 19, 11)
gmsh.model.occ.addLine(4, 20, 12)
gmsh.model.occ.addLine(9, 2, 13)
gmsh.model.occ.addLine(10, 3, 14)
gmsh.model.occ.addLine(11, 4, 15)
gmsh.model.occ.addLine(12, 1, 16)
gmsh.model.occ.addLine(13, 6, 17)
gmsh.model.occ.addLine(14, 7, 18)
gmsh.model.occ.addLine(15, 8, 19)
gmsh.model.occ.addLine(16, 5, 20)
gmsh.model.occ.addLine(17, 6, 21)
gmsh.model.occ.addLine(18, 5, 22)
gmsh.model.occ.addLine(19, 7, 23)
gmsh.model.occ.addLine(20, 8, 24)

gmsh.model.occ.addCurveLoop([1, 2, 3, 4], 1)
gmsh.model.occ.addCurveLoop([5, 6, 7, 8], 2)
gmsh.model.occ.addCurveLoop([10, 5, 9, 1], 3)
gmsh.model.occ.addCurveLoop([12, 7, 11, 3], 4)
gmsh.model.occ.addCurveLoop([9, 6, 11, 2], 5)
gmsh.model.occ.addCurveLoop([10, 8, 12, 4], 6)

gmsh.model.occ.addSurfaceFilling(1, 1)
gmsh.model.occ.addSurfaceFilling(2, 2)
gmsh.model.occ.addSurfaceFilling(3, 3)
gmsh.model.occ.addSurfaceFilling(4, 4)
gmsh.model.occ.addSurfaceFilling(5, 5)
gmsh.model.occ.addSurfaceFilling(6, 6)

gmsh.model.occ.addLine(9, 15, 25)
gmsh.model.occ.addLine(10, 16, 26)
gmsh.model.occ.addLine(11, 13, 27)
gmsh.model.occ.addLine(12, 14, 28)
gmsh.model.occ.addLine(17, 20, 29)
gmsh.model.occ.addLine(18, 19, 30)
gmsh.model.occ.addSpline([1,9,2],1)
gmsh.model.occ.addSpline([2,10,3],2)
gmsh.model.occ.addSpline([3,11,4],3)
gmsh.model.occ.addSpline([4,12,1],4)
gmsh.model.occ.addSpline([5,13,6],5)
gmsh.model.occ.addSpline([6,14,7],6)
gmsh.model.occ.addSpline([7,15,8],7)
gmsh.model.occ.addSpline([8,16,5],8)
gmsh.model.occ.addSpline([2,17,6],9)
gmsh.model.occ.addSpline([1,18,5],10)
gmsh.model.occ.addSpline([3,19,7],11)
gmsh.model.occ.addSpline([4,20,8],12)


gmsh.model.occ.addCurveLoop([1,2,3,4],1)
gmsh.model.occ.addCurveLoop([5,6,7,8],2)
gmsh.model.occ.addCurveLoop([10,5,9,1],3)
gmsh.model.occ.addCurveLoop([12,7,11,3],4)
gmsh.model.occ.addCurveLoop([9,6,11,2],5)
gmsh.model.occ.addCurveLoop([10,8,12,4],6)

gmsh.model.occ.addSurfaceFilling(1,1)
gmsh.model.occ.addSurfaceFilling(2,2)
gmsh.model.occ.addSurfaceFilling(3,3)
gmsh.model.occ.addSurfaceFilling(4,4)
gmsh.model.occ.addSurfaceFilling(5,5)
gmsh.model.occ.addSurfaceFilling(6,6)

gmsh.model.occ.addLine(9,15,13)
gmsh.model.occ.addLine(10,16,14)
gmsh.model.occ.addLine(11,13,15)
gmsh.model.occ.addLine(12,14,16)
gmsh.model.occ.addLine(17,20,17)
gmsh.model.occ.addLine(18,19,18)
gmsh.model.occ.synchronize()
"""
def assemble_models(model_name="Assembly",model_names=[],lc=None,gui=False):
    if not gmsh.isInitialized():
        gmsh.initialize()
    gmsh.model.add(model_name)
    for model_name in model_names:
        gmsh.merge(model_name)
    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()

    # create map containing information on mesh periodicity and physical groups
    group_map = {"POINT_1": [1],
                 "POINT_2": [2],
                 "POINT_3": [3],
                 "POINT_4": [4],
                 "POINT_5": [5],
                 "POINT_6": [6],
                 "POINT_7": [7],
                 "POINT_8": [8],
                 "EDGE_1": [1],
                 "EDGE_2": [2],
                 "EDGE_3": [3],
                 "EDGE_4": [4],
                 "EDGE_5": [5],
                 "EDGE_6": [6],
                 "EDGE_7": [7],
                 "EDGE_8": [8],
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
    gmsh.model.mesh.reclassifyNodes()

    if gui:
        gmsh.fltk.run()

    return model_name,group_map

#assemble_models(model_names=["C:/Users/lemji/OneDrive/Bureau/khedma/BCC","C:/Users/lemji/OneDrive/Bureau/khedma/FCC"],lc=5,gui=True)
