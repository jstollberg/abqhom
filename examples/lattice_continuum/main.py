import numpy as np
import gmsh


def create_geometry(model_name, dx, dy, dz, radius_cube, radius_bcc):

    # Initialize Gmsh and add a new model
    if not gmsh.is_initialized():
        gmsh.initialize()
    gmsh.model.add(model_name)

    occ = gmsh.model.occ
    gmsh.option.set_number("Geometry.Tolerance", 1e-5)

    # Construct the outer cube
    cube_strut_1 = occ.add_cylinder(-0.5*dx, -0.5*dy, -0.5*dz, dx, 0.0, 0.0,
                                    radius_cube)
    cube_strut_2 = occ.add_cylinder(0.5*dx, -0.5*dy, -0.5*dz, 0.0, dy, 0.0,
                                    radius_cube)
    cube_strut_3 = occ.add_cylinder(0.5*dx, 0.5*dy, -0.5*dz, -dx, 0.0, 0.0,
                                    radius_cube)
    cube_strut_4 = occ.add_cylinder(-0.5*dx, 0.5*dy, -0.5*dz, 0.0, -dy, 0.0,
                                    radius_cube)
    cube_strut_5 = occ.add_cylinder(-0.5*dx, -0.5*dy, 0.5*dz, dx, 0.0, 0.0,
                                    radius_cube)
    cube_strut_6 = occ.add_cylinder(0.5*dx, -0.5*dy, 0.5*dz, 0.0, dy, 0.0,
                                    radius_cube)
    cube_strut_7 = occ.add_cylinder(0.5*dx, 0.5*dy, 0.5*dz, -dx, 0.0, 0.0,
                                    radius_cube)
    cube_strut_8 = occ.add_cylinder(-0.5*dx, 0.5*dy, 0.5*dz, 0.0, -dy, 0.0,
                                    radius_cube)
    cube_strut_9 = occ.add_cylinder(-0.5*dx, -0.5*dy, -0.5*dz, 0.0, 0.0, dz,
                                    radius_cube)
    cube_strut_10 = occ.add_cylinder(0.5*dx, -0.5*dy, -0.5*dz, 0.0, 0.0, dz,
                                     radius_cube)
    cube_strut_11 = occ.add_cylinder(0.5*dx, 0.5*dy, -0.5*dz, 0.0, 0.0, dz,
                                     radius_cube)
    cube_strut_12 = occ.add_cylinder(-0.5*dx, 0.5*dy, -0.5*dz, 0.0, 0.0, dz,
                                     radius_cube)

    cube_dimtags = []
    for tag in [cube_strut_1, cube_strut_2, cube_strut_3, cube_strut_4,
                cube_strut_5, cube_strut_6, cube_strut_7, cube_strut_8,
                cube_strut_9, cube_strut_10, cube_strut_11, cube_strut_12]:
        cube_dimtags.append((3, tag))
    cube_struts, _ = occ.fuse(cube_dimtags, cube_dimtags)

    # Construct the BCC struts
    bcc_strut_1 = occ.add_cylinder(-0.5*dx, -0.5*dy, -0.5*dz, dx, dy, dz,
                                   radius_bcc)
    bcc_strut_2 = occ.add_cylinder(0.5*dx, -0.5*dy, -0.5*dz, -dx, dy, dz,
                                   radius_bcc)
    bcc_strut_3 = occ.add_cylinder(0.5*dx, 0.5*dy, -0.5*dz, -dx, -dy, dz,
                                   radius_bcc)
    bcc_strut_4 = occ.add_cylinder(-0.5*dx, 0.5*dy, -0.5*dz, dx, -dy, dz,
                                   radius_bcc)

    bcc_dimtags = []
    for tag in [bcc_strut_1, bcc_strut_2, bcc_strut_3, bcc_strut_4]:
        bcc_dimtags.append((3, tag))
    bcc_struts, _ = occ.fuse(bcc_dimtags, bcc_dimtags)

    # Perform boolean operations to cut the struts to a cube
    bounding_box = occ.add_box(-0.5*dx, -0.5*dy, -0.5*dz, dx, dy, dz)
    cube_struts, _ = occ.intersect(cube_struts, [(3, bounding_box)],
                                   removeTool=False)
    bcc_struts, _ = occ.intersect(bcc_struts, [(3, bounding_box)])
    bcc_struts, _ = occ.cut(bcc_struts, cube_struts, removeTool=False)
    gmsh.model.occ.synchronize()

    # Add physical groups
    cube = gmsh.model.add_physical_group(3, [cube_struts[0][1]])
    gmsh.model.set_physical_name(3, cube, "CUBIC")
    bcc = gmsh.model.add_physical_group(3, [bcc_struts[0][1]])
    gmsh.model.set_physical_name(3, bcc, "BCC")



from abqhom.gmshmodel import Boundary, GmshRVE



dx = dy = dz = 10.0
radius_cube = 1.0
radius_bcc = 2.0
model_name = "RUC"

gmsh.initialize()


create_geometry(model_name, dx, dy, dz, radius_cube, radius_bcc)
unit_cell = GmshRVE(model_name)
unit_cell.generate_mesh(element_size=0.5)
# unit_cell.build_element_sets()
# unit_cell._periodic=False

# #

# unit_cell.build_node_sets()
# print(unit_cell.node_sets["FACE_NODES"])
# print(unit_cell.node_sets["EDGE_4"])
unit_cell.show()






# gmsh.finalize()
