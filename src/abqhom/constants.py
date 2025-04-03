"""
Collection of constants for the abqhom package.

"""

# Physical groups on the RVE boundary
POINT_GROUPS = ["POINT_1", "POINT_2", "POINT_3", "POINT_4",
                "POINT_5", "POINT_6", "POINT_7", "POINT_8"]
EDGE_GROUPS = ["EDGE_1", "EDGE_2", "EDGE_3", "EDGE_4", "EDGE_5", "EDGE_6",
               "EDGE_7", "EDGE_8", "EDGE_9", "EDGE_10", "EDGE_11", "EDGE_12"]
FACE_GROUPS = ["FACE_1", "FACE_2", "FACE_3", "FACE_4", "FACE_5", "FACE_6"]

# Periodic pairs of groups, keys are masters, values are slaves
PERIODIC_PAIRS = {"EDGE_1":  "EDGE_3",
                  "EDGE_3":  "EDGE_7",
                  "EDGE_7":  "EDGE_5",
                  "EDGE_2":  "EDGE_4",
                  "EDGE_4":  "EDGE_8",
                  "EDGE_8":  "EDGE_6",
                  "EDGE_9":  "EDGE_10",
                  "EDGE_10": "EDGE_11",
                  "EDGE_11": "EDGE_12",
                  "FACE_1":  "FACE_2",
                  "FACE_3":  "FACE_4",
                  "FACE_5":  "FACE_6"}

# Groups forming the boundary of a higher-dimensional group
BOUNDARY_CONNECTIVITY = {"POINT_1": [],
                         "POINT_2": [],
                         "POINT_3": [],
                         "POINT_4": [],
                         "POINT_5": [],
                         "POINT_6": [],
                         "POINT_7": [],
                         "POINT_8": [],
                         "EDGE_1":  ["POINT_1", "POINT_2"],
                         "EDGE_2":  ["POINT_2", "POINT_3"],
                         "EDGE_3":  ["POINT_3", "POINT_4"],
                         "EDGE_4":  ["POINT_4", "POINT_1"],
                         "EDGE_5":  ["POINT_5", "POINT_6"],
                         "EDGE_6":  ["POINT_6", "POINT_7"],
                         "EDGE_7":  ["POINT_7", "POINT_8"],
                         "EDGE_8":  ["POINT_8", "POINT_5"],
                         "EDGE_9":  ["POINT_1", "POINT_5"],
                         "EDGE_10": ["POINT_2", "POINT_6"],
                         "EDGE_11": ["POINT_3", "POINT_7"],
                         "EDGE_12": ["POINT_4", "POINT_8"],
                         "FACE_1":  ["EDGE_1", "EDGE_2", "EDGE_3", "EDGE_4"],
                         "FACE_2":  ["EDGE_5", "EDGE_6", "EDGE_7", "EDGE_8"],
                         "FACE_3":  ["EDGE_1", "EDGE_9", "EDGE_5", "EDGE_10"],
                         "FACE_4":  ["EDGE_3", "EDGE_12", "EDGE_7", "EDGE_11"],
                         "FACE_5":  ["EDGE_9", "EDGE_4", "EDGE_12", "EDGE_8"],
                         "FACE_6":  ["EDGE_10", "EDGE_2", "EDGE_11", "EDGE_6"]}

# Default map from Gmsh element types to Abaqus element types
ABAQUS_ELEMENT_TYPES = {1:  "B23",    # 2-node line
                        2:  "CPE3",   # 3-node triangle
                        3:  "CPE4",   # 4-node quadrangle
                        4:  "C3D4",   # 4-node tetrahedron
                        5:  "C3D8",   # 8-node hexahedron
                        6:  "C3D6",   # 6-node prism
                        7:  "C3D5",   # 5-node pyramid
                        8:  "B22",    # 3-node second-order line
                        9:  "CPE6",   # 6-node second order triangle
                        11: "C3D10",  # 10-node second-order tetrahedron
                        16: "CPE8",   # 8-node second-order quadrangle
                        17: "C3D20",  # 20-node second-order hexahedron
                        18: "C3D15",  # 15-node second-order prism
                        }

# Map from the order of the connectivity list from Gmsh to Abaqus
ABAQUS_CONNECTIVITY_ORDER = {
    1:  [0, 1],                               # 2-node line
    2:  [0, 1, 2],                            # 3-node triangle
    3:  [0, 1, 2, 3],                         # 4-node quadrangle
    4:  [0, 1, 2, 3],                         # 4-node tetrahedron
    5:  [0, 1, 2, 3, 4, 5, 6, 7],             # 8-node hexahedron
    6:  [0, 1, 2, 3, 4, 5],                   # 6-node prism
    7:  [0, 1, 2, 3, 4],                      # 5-node pyramid
    8:  [0, 2, 1],                            # 3-node 2nd line
    9:  [0, 1, 2, 3, 4, 5],                   # 6-node 2nd triangle
    11: [0, 1, 2, 3, 4, 5, 6, 7, 9, 8],       # 10-node 2nd tetrahedron
    16: [0, 1, 2, 3, 4, 5, 6, 7],             # 8-node 2nd quadrangle
    17: [0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 13,
         9, 16, 18, 19, 17, 10, 12, 14, 15],  # 20-node 2nd hexahedron
    18: [0, 1, 2, 3, 4, 5, 7, 9,
         7, 12, 14, 12, 8, 10, 11],           # 15-node 2nd prism
                             }
