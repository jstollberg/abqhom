# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-only

"""
Collection of constants for the abqhom package.

"""

try:
    import paraqus  # noqa: F401
    PARAQUS_INSTALLED = True
except ImportError:
    PARAQUS_INSTALLED = False

# Physical groups on the RVE boundary
POINT_GROUPS = ["POINT_1", "POINT_2", "POINT_3", "POINT_4",
                "POINT_5", "POINT_6", "POINT_7", "POINT_8"]
EDGE_GROUPS = ["EDGE_1", "EDGE_2", "EDGE_3", "EDGE_4", "EDGE_5", "EDGE_6",
               "EDGE_7", "EDGE_8", "EDGE_9", "EDGE_10", "EDGE_11", "EDGE_12"]
FACE_GROUPS = ["FACE_1", "FACE_2", "FACE_3", "FACE_4", "FACE_5", "FACE_6"]

# Periodic pairs of groups, with secondary groups as keys and primary groups
# as values
PERIODIC_PAIRS = {"POINT_2": "POINT_1",
                  "POINT_3": "POINT_1",
                  "POINT_4": "POINT_1",
                  "POINT_5": "POINT_1",
                  "POINT_6": "POINT_1",
                  "POINT_7": "POINT_1",
                  "POINT_8": "POINT_1",
                  "EDGE_3":  "EDGE_1",
                  "EDGE_5":  "EDGE_1",
                  "EDGE_7":  "EDGE_1",
                  "EDGE_2":  "EDGE_4",
                  "EDGE_6":  "EDGE_4",
                  "EDGE_8":  "EDGE_4",
                  "EDGE_10": "EDGE_9",
                  "EDGE_11": "EDGE_9",
                  "EDGE_12": "EDGE_9",
                  "FACE_2":  "FACE_1",
                  "FACE_4":  "FACE_3",
                  "FACE_6":  "FACE_5"}

# Translation axes from each primary group to its secondary partner
SECONDARY_GROUPS_OFFSET_AXES = {"POINT_2": [0],
                                "POINT_3": [0, 1],
                                "POINT_4": [1],
                                "POINT_5": [2],
                                "POINT_6": [0, 2],
                                "POINT_7": [0, 1, 2],
                                "POINT_8": [1, 2],
                                "EDGE_3":  [1],
                                "EDGE_5":  [2],
                                "EDGE_7":  [1, 2],
                                "EDGE_2":  [0],
                                "EDGE_6":  [0, 2],
                                "EDGE_8":  [2],
                                "EDGE_10": [0],
                                "EDGE_11": [0, 1],
                                "EDGE_12": [1],
                                "FACE_2":  [2],
                                "FACE_4":  [1],
                                "FACE_6":  [0]}

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

# Default map from Gmsh element types to Abaqus element types for 2D models
ABAQUS_ELEMENT_TYPES_2D = {1:  "B23",    # 2-node line
                           2:  "CPE3",   # 3-node triangle
                           3:  "CPE4",   # 4-node quadrangle
                           8:  "B22",    # 3-node second-order line
                           9:  "CPE6",   # 6-node second-order triangle
                           16: "CPE8",   # 8-node second-order quadrangle
                           }

# Default map from Gmsh element types to Abaqus element types for 3D models
ABAQUS_ELEMENT_TYPES_3D = {1:  "B33",     # 2-node line
                           2:  "S3",      # 3-node triangle
                           3:  "S4",      # 4-node quadrangle
                           4:  "C3D4",    # 4-node tetrahedron
                           5:  "C3D8",    # 8-node hexahedron
                           6:  "C3D6",    # 6-node prism
                           7:  "C3D5",    # 5-node pyramid
                           8:  "B32",     # 3-node second-order line
                           9:  "STRI65",  # 6-node second-order triangle
                           11: "C3D10",   # 10-node second-order tetrahedron
                           16: "S8R5",    # 8-node second-order quadrangle
                           17: "C3D20",   # 20-node second-order hexahedron
                           18: "C3D15",   # 15-node second-order prism
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
    18: [0, 1, 2, 3, 4, 5, 6, 9,
         7, 12, 14, 13, 8, 10, 11],           # 15-node 2nd prism
                             }

# Abaqus element types with rotational degrees of freedom
ELEMENTS_ROTATIONAL_DOFS = {
    # Beams in a plane
    "B21", "B21H", "B22", "B22H", "B23", "B23H", "PIPE21", "PIPE21H", "PIPE22",
    "PIPE22H",

    # Beams in space
    "B31", "B31H", "B32", "B32H", "B33", "B33H", "PIPE31", "PIPE31H", "PIPE32",
    "PIPE32H",

    # Frame elements
    "FRAME2D", "FRAME3D",

    # Elbow elements
    "ELBOW31", "ELBOW32", "ELBOW31B", "ELBOW31C",

    # Shell elements
    "STRI3", "S3", "S3R", "S3RS", "STRI65", "S4", "S4R", "S4RS", "S4RSW",
    "S4R5", "S8R", "S8R5", "S9R5"
    }

# Abaqus element types with warping degrees of freedom
ELEMENTS_WARPING_DOFS = {
    "B31OS", "B31OSH", "B32OS", "B32OSH"
    }

# Abaqus element types that are not supported by abqhom
ELEMENTS_NOT_SUPPORTED = {
    "STRI65", "S4R5", "S8R5", "S9R5"
    }
