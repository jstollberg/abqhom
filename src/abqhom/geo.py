import gmsh
from abqhom.RVE import (add_boundary_groups, add_periodic_constraints,
                        create_boundary_sets)
from abqhom.examples import *

"""
trying to work on the assembly method here
"""

def BBC_geo(model_name="BCC", x=0.0, y=0.0, z=0.0, dx=1.0, dy=1.0,
                   dz=1.0, lc=None, gui=False):


    # create Geometry

    gmsh.model.occ.addLine(4,6)
    gmsh.model.occ.addLine(3,5)
    gmsh.model.occ.addLine(1,7)
    gmsh.model.occ.addLine(2,8)

def FCC_geo(model_name="FCC", x=0.0, y=0.0, z=0.0, dx=1.0, dy=1.0,
                   dz=1.0, lc=None, gui=False):

    Dx, Dy, Dz = dx / 2, dy / 2, dz / 2

    tag_p9=gmsh.model.occ.addPoint(x+Dx,y,z,lc)
    tag_p10=gmsh.model.occ.addPoint(x-Dx,y,z,lc)
    tag_p11=gmsh.model.occ.addPoint(x,y+Dy,z,lc)
    tag_p12=gmsh.model.occ.addPoint(x,y-Dy,z,lc)
    tag_p13=gmsh.model.occ.addPoint(x,y,z+Dz,lc)
    tag_p14=gmsh.model.occ.addPoint(x,y,z-Dz,lc)

    tag_l13=gmsh.model.occ.addLine(tag_p9,tag_p10)
    tag_l14=gmsh.model.occ.addLine(tag_p11,tag_p12)
    tag_l15=gmsh.model.occ.addLine(tag_p13,tag_p14)


    gmsh.model.addPhysicalGroup(0, [tag_p9])
    gmsh.model.addPhysicalGroup(0, [tag_p10])
    gmsh.model.addPhysicalGroup(0, [tag_p11])
    gmsh.model.addPhysicalGroup(0, [tag_p12])
    gmsh.model.addPhysicalGroup(0, [tag_p13])
    gmsh.model.addPhysicalGroup(0, [tag_p14])

    gmsh.model.addPhysicalGroup(1, [tag_l13])
    gmsh.model.addPhysicalGroup(1, [tag_l14])
    gmsh.model.addPhysicalGroup(1, [tag_l15])
def AFFC_geo(model_name="AFCC", x=0.0, y=0.0, z=0.0, dx=1.0, dy=1.0,
                   dz=1.0, lc=None, gui=False):

    gmsh.model.occ.addLine(1,6)
    gmsh.model.occ.addLine(1,3)
    gmsh.model.occ.addLine(1,8)
    gmsh.model.occ.addLine(2,4)
    gmsh.model.occ.addLine(2,5)
    gmsh.model.occ.addLine(2,7)
    gmsh.model.occ.addLine(3, 6)
    gmsh.model.occ.addLine(3, 8)
    gmsh.model.occ.addLine(4,5)
    gmsh.model.occ.addLine(4,7)
    gmsh.model.occ.addLine(5,7)
    gmsh.model.occ.addLine(6,8)

def diamond_geo(model_name="Diamond", x=0.0, y=0.0, z=0.0, dx=1.0, dy=1.0,
                   dz=1.0, lc=None, gui=False):

    Dx, Dy, Dz = dx / 2, dy / 2, dz / 2

    tag_p9=gmsh.model.occ.addPoint(x + dx / 2, y, z,lc)
    tag_p10=gmsh.model.occ.addPoint(x - dx / 2, y, z,lc)
    tag_p11=gmsh.model.occ.addPoint(x, y + dy / 2, z,lc)
    tag_p12=gmsh.model.occ.addPoint(x, y - dy / 2, z,lc)
    tag_p13=gmsh.model.occ.addPoint(x, y, z + dz / 2,lc)
    tag_p14=gmsh.model.occ.addPoint(x, y, z - dz / 2,lc)

    tag_l13=gmsh.model.occ.addLine(9,11)
    tag_l14=gmsh.model.occ.addLine(9,12)
    tag_l15=gmsh.model.occ.addLine(10,11)
    tag_l16=gmsh.model.occ.addLine(10,12)
    tag_l17=gmsh.model.occ.addLine(11,13)
    tag_l18=gmsh.model.occ.addLine(11,14)
    tag_l19=gmsh.model.occ.addLine(12,13)
    tag_l20=gmsh.model.occ.addLine(12,14)
    tag_l21=gmsh.model.occ.addLine(11,12)




    gmsh.model.addPhysicalGroup(0, [tag_p9])
    gmsh.model.addPhysicalGroup(0, [tag_p10])
    gmsh.model.addPhysicalGroup(0, [tag_p11])
    gmsh.model.addPhysicalGroup(0, [tag_p12])
    gmsh.model.addPhysicalGroup(0, [tag_p13])
    gmsh.model.addPhysicalGroup(0, [tag_p14])

    gmsh.model.addPhysicalGroup(1, [tag_l13])
    gmsh.model.addPhysicalGroup(1, [tag_l14])
    gmsh.model.addPhysicalGroup(1, [tag_l15])
    gmsh.model.addPhysicalGroup(1, [tag_l16])
    gmsh.model.addPhysicalGroup(1, [tag_l17])
    gmsh.model.addPhysicalGroup(1, [tag_l18])
    gmsh.model.addPhysicalGroup(1, [tag_l19])
    gmsh.model.addPhysicalGroup(1, [tag_l20])
    gmsh.model.addPhysicalGroup(1, [tag_l21])
def AAC_geo(model_name="AAC", x=0.0, y=0.0, z=0.0, dx=1.0, dy=1.0,
                   dz=1.0, lc=None, gui=False):
    Dx, Dy, Dz = dx / 2, dy / 2, dz / 2

    tag_p9=gmsh.model.occ.addPoint(x - Dx, y - Dy, z, lc,)
    tag_p10=gmsh.model.occ.addPoint(x - Dx, y, z - Dz, lc)
    tag_p11=gmsh.model.occ.addPoint(x - Dx, y + Dy, z, lc)
    tag_p12=gmsh.model.occ.addPoint(x - Dx, y, z + Dz, lc)
    tag_p13=gmsh.model.occ.addPoint(x + Dx, y - Dy, z, lc)
    tag_p14=gmsh.model.occ.addPoint(x + Dx, y, z - Dz, lc)
    tag_p15=gmsh.model.occ.addPoint(x + Dx, y + Dy, z, lc)
    tag_p16=gmsh.model.occ.addPoint(x + Dx, y, z + Dz, lc)
    tag_p17=gmsh.model.occ.addPoint(x, y - Dy, z - Dz, lc)
    tag_p18=gmsh.model.occ.addPoint(x, y - Dy, z + Dz, lc)
    tag_p19=gmsh.model.occ.addPoint(x, y + Dy, z - Dz, lc)
    tag_p20=gmsh.model.occ.addPoint(x, y + Dy, z + Dz, lc)


    gmsh.model.occ.addSpline([1, tag_p9, 2],1)
    gmsh.model.occ.addSpline([2, tag_p10, 3],2)
    gmsh.model.occ.addSpline([3, tag_p11, 4],3)
    gmsh.model.occ.addSpline([4, tag_p12, 1],4)
    gmsh.model.occ.addSpline([5, tag_p13, 6],5)
    gmsh.model.occ.addSpline([6, tag_p14, 7],6)
    gmsh.model.occ.addSpline([7, tag_p15, 8],7)
    gmsh.model.occ.addSpline([8, tag_p16, 5],8)
    gmsh.model.occ.addSpline([2, tag_p17, 6],9)
    gmsh.model.occ.addSpline([1,tag_p18, 5],10)
    gmsh.model.occ.addSpline([3, tag_p19, 7],11)
    gmsh.model.occ.addSpline([4, tag_p20, 8],12)

    tag_l13=gmsh.model.occ.addLine(tag_p9, tag_p15)
    tag_l14=gmsh.model.occ.addLine(tag_p10, tag_p16)
    tag_l15=gmsh.model.occ.addLine(tag_p11, tag_p13)
    tag_l16=gmsh.model.occ.addLine(tag_p12, tag_p14)
    tag_l17=gmsh.model.occ.addLine(tag_p17, tag_p20)
    tag_l18=gmsh.model.occ.addLine(tag_p18, tag_p19)



    gmsh.model.occ.fragment([(1, 1)], [(1, tag_l13)])
    gmsh.model.occ.fragment([(1, 2)], [(1, tag_l14)])
    gmsh.model.occ.fragment([(1, 3)], [(1, tag_l15)])
    gmsh.model.occ.fragment([(1, 4)], [(1, tag_l16)])
    gmsh.model.occ.fragment([(1, 5)], [(1, tag_l15)])
    gmsh.model.occ.fragment([(1, 6)], [(1, tag_l16)])
    gmsh.model.occ.fragment([(1, 7)], [(1, tag_l13)])
    gmsh.model.occ.fragment([(1, 8)], [(1, tag_l14)])
    gmsh.model.occ.fragment([(1, 9)], [(1, tag_l17)])
    gmsh.model.occ.fragment([(1, 10)], [(1, tag_l18)])
    gmsh.model.occ.fragment([(1, 11)], [(1, tag_l18)])
    gmsh.model.occ.fragment([(1, 12)], [(1, tag_l17)])

    gmsh.model.addPhysicalGroup(1, [tag_l13])
    gmsh.model.addPhysicalGroup(1, [tag_l14])
    gmsh.model.addPhysicalGroup(1, [tag_l15])
    gmsh.model.addPhysicalGroup(1, [tag_l16])
    gmsh.model.addPhysicalGroup(1, [tag_l17])
    gmsh.model.addPhysicalGroup(1, [tag_l18])

    gmsh.model.addPhysicalGroup(0, [tag_p9])
    gmsh.model.addPhysicalGroup(0, [tag_p10])
    gmsh.model.addPhysicalGroup(0, [tag_p11])
    gmsh.model.addPhysicalGroup(0, [tag_p12])
    gmsh.model.addPhysicalGroup(0, [tag_p13])
    gmsh.model.addPhysicalGroup(0, [tag_p14])
    gmsh.model.addPhysicalGroup(0, [tag_p15])
    gmsh.model.addPhysicalGroup(0, [tag_p16])
    gmsh.model.addPhysicalGroup(0, [tag_p17])
    gmsh.model.addPhysicalGroup(0, [tag_p18])
    gmsh.model.addPhysicalGroup(0, [tag_p19])
    gmsh.model.addPhysicalGroup(0, [tag_p20])

def Central_cube_geo(model_name="CCube", x=0.0, y=0.0, z=0.0, dx=1.0, dy=1.0,
                   dz=1.0,ratio_x=1/2,ratio_y=1/2,ratio_z=1/2, lc=None, gui=False):
    Dx, Dy, Dz = dx*ratio_x / 2, dy*ratio_y / 2, dz*ratio_z / 2

    tag_p1=gmsh.model.occ.addPoint(x - Dx, y - Dy, z + Dz, lc)
    tag_p2=gmsh.model.occ.addPoint(x - Dx, y - Dy, z - Dz, lc)
    tag_p3=gmsh.model.occ.addPoint(x - Dx, y + Dy, z - Dz, lc)
    tag_p4=gmsh.model.occ.addPoint(x - Dx, y + Dy, z + Dz, lc)
    tag_p5=gmsh.model.occ.addPoint(x + Dx, y - Dy, z + Dz, lc)
    tag_p6=gmsh.model.occ.addPoint(x + Dx, y - Dy, z - Dz, lc)
    tag_p7=gmsh.model.occ.addPoint(x + Dx, y + Dy, z - Dz, lc)
    tag_p8=gmsh.model.occ.addPoint(x + Dx, y + Dy, z + Dz, lc)

    tag_l1=gmsh.model.occ.addLine(tag_p1, tag_p2)
    tag_l2=gmsh.model.occ.addLine(tag_p2, tag_p3)
    tag_l3=gmsh.model.occ.addLine(tag_p3, tag_p4)
    tag_l4=gmsh.model.occ.addLine(tag_p4, tag_p1)
    tag_l5=gmsh.model.occ.addLine(tag_p5, tag_p6)
    tag_l6=gmsh.model.occ.addLine(tag_p6, tag_p7)
    tag_l7=gmsh.model.occ.addLine(tag_p7, tag_p8)
    tag_l8=gmsh.model.occ.addLine(tag_p8, tag_p5)
    tag_l9=gmsh.model.occ.addLine(tag_p2, tag_p6)
    tag_l10=gmsh.model.occ.addLine(tag_p1, tag_p5)
    tag_l11=gmsh.model.occ.addLine(tag_p3, tag_p7)
    tag_l12=gmsh.model.occ.addLine(tag_p4, tag_p8)

    gmsh.model.occ.addLine(1,tag_p1)
    gmsh.model.occ.addLine(2,tag_p2)
    gmsh.model.occ.addLine(3,tag_p3)
    gmsh.model.occ.addLine(4,tag_p4)
    gmsh.model.occ.addLine(5,tag_p5)
    gmsh.model.occ.addLine(6,tag_p6)
    gmsh.model.occ.addLine(7,tag_p7)
    gmsh.model.occ.addLine(8,tag_p8)

    gmsh.model.addPhysicalGroup(1, [tag_l1])
    gmsh.model.addPhysicalGroup(1, [tag_l2])
    gmsh.model.addPhysicalGroup(1, [tag_l3])
    gmsh.model.addPhysicalGroup(1, [tag_l4])
    gmsh.model.addPhysicalGroup(1, [tag_l5])
    gmsh.model.addPhysicalGroup(1, [tag_l6])
    gmsh.model.addPhysicalGroup(1, [tag_l7])
    gmsh.model.addPhysicalGroup(1, [tag_l8])
    gmsh.model.addPhysicalGroup(1, [tag_l9])
    gmsh.model.addPhysicalGroup(1, [tag_l10])
    gmsh.model.addPhysicalGroup(1, [tag_l11])
    gmsh.model.addPhysicalGroup(1, [tag_l12])

    gmsh.model.addPhysicalGroup(0, [tag_p1])
    gmsh.model.addPhysicalGroup(0, [tag_p2])
    gmsh.model.addPhysicalGroup(0, [tag_p3])
    gmsh.model.addPhysicalGroup(0, [tag_p4])
    gmsh.model.addPhysicalGroup(0, [tag_p5])
    gmsh.model.addPhysicalGroup(0, [tag_p6])
    gmsh.model.addPhysicalGroup(0, [tag_p7])
    gmsh.model.addPhysicalGroup(0, [tag_p8])




def assemble(model_name="Assembly", x=0.0, y=0.0, z=0.0, dx=1.0, dy=1.0,
             dz=1.0, lc=None, gui=False, parts=None ,ACC=False):
    if parts is None:
        parts = []
    if ACC:
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

        gmsh.model.occ.addPoint(x - Dx, y - Dy, z, lc, 9)
        gmsh.model.occ.addPoint(x - Dx, y, z - Dz, lc, 10)
        gmsh.model.occ.addPoint(x - Dx, y + Dy, z, lc, 11)
        gmsh.model.occ.addPoint(x - Dx, y, z + Dz, lc, 12)
        gmsh.model.occ.addPoint(x + Dx, y - Dy, z, lc, 13)
        gmsh.model.occ.addPoint(x + Dx, y, z - Dz, lc, 14)
        gmsh.model.occ.addPoint(x + Dx, y + Dy, z, lc, 15)
        gmsh.model.occ.addPoint(x + Dx, y, z + Dz, lc, 16)
        gmsh.model.occ.addPoint(x, y - Dy, z - Dz, lc, 17)
        gmsh.model.occ.addPoint(x, y - Dy, z + Dz, lc, 18)
        gmsh.model.occ.addPoint(x, y + Dy, z - Dz, lc, 19)
        gmsh.model.occ.addPoint(x, y + Dy, z + Dz, lc, 20)

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

        for p in parts:
            if p=="BCC":
                BBC_geo(x=x, y=y, z=z, dx=dx, dy=dy, dz=dz,lc=lc)
            elif p=="FCC":
                FCC_geo(x=x,y=y,z=z, dx=dx, dy=dy, dz=dz,lc=lc)
            elif p=="AFCC":
                AFFC_geo(x=x,y=y,z=z, dx=dx, dy=dy, dz=dz,lc=lc)
            elif p=="Diamond":
                diamond_geo(x=x,y=y,z=z, dx=dx, dy=dy, dz=dz,lc=lc)
            elif p=="Cubic":
                Central_cube_geo(x=x,y=y,z=z, dx=dx, dy=dy, dz=dz,lc=lc)

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

        gmsh.model.occ.addSurfaceLoop([1, 2, 3, 4, 5, 6], 1)
        gmsh.model.occ.addVolume([1], 1)

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

        gmsh.model.addPhysicalGroup(1, [13])
        gmsh.model.addPhysicalGroup(1, [14])
        gmsh.model.addPhysicalGroup(1, [15])
        gmsh.model.addPhysicalGroup(1, [16])
        gmsh.model.addPhysicalGroup(1, [17])
        gmsh.model.addPhysicalGroup(1, [18])

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
    else:
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

        for p in parts:
            if p=="BCC":
                BBC_geo(x=x, y=y, z=z, dx=dx, dy=dy, dz=dz,lc=lc)
            elif p=="FCC":
                FCC_geo(x=x,y=y,z=z, dx=dx, dy=dy, dz=dz,lc=lc)
            elif p=="AFCC":
                AFFC_geo(x=x,y=y,z=z, dx=dx, dy=dy, dz=dz,lc=lc)
            elif p=="Diamond":
                diamond_geo(x=x,y=y,z=z, dx=dx, dy=dy, dz=dz,lc=lc)
            elif p=="Cubic":
                Central_cube_geo(x=x,y=y,z=z, dx=dx, dy=dy, dz=dz,lc=lc)


        Boundary=[1,2,3,4,5,6,7,8,9,10,11,12]




        gmsh.model.occ.addCurveLoop([Boundary[0], Boundary[1], Boundary[2], Boundary[3]], 1)
        gmsh.model.occ.addCurveLoop([Boundary[4], Boundary[5], Boundary[6], Boundary[7]], 2)
        gmsh.model.occ.addCurveLoop([Boundary[9], Boundary[4], Boundary[8], Boundary[0]], 3)
        gmsh.model.occ.addCurveLoop([Boundary[11], Boundary[6], Boundary[10], Boundary[2]], 4)
        gmsh.model.occ.addCurveLoop([Boundary[8], Boundary[5], Boundary[10], Boundary[1]], 5)
        gmsh.model.occ.addCurveLoop([Boundary[9], Boundary[7], Boundary[11], Boundary[3]], 6)

        gmsh.model.occ.addSurfaceFilling(1, 1)
        gmsh.model.occ.addSurfaceFilling(2, 2)
        gmsh.model.occ.addSurfaceFilling(3, 3)
        gmsh.model.occ.addSurfaceFilling(4, 4)
        gmsh.model.occ.addSurfaceFilling(5, 5)
        gmsh.model.occ.addSurfaceFilling(6, 6)

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


        if gui:
            gmsh.fltk.run()

        return model_name, group_map
assemble(lc=5,gui=True,parts=["BCC","FCC"],ACC=True)
