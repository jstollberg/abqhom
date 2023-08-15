import numpy as np
import gmsh
import os

GROUP_TO_DIM = {"POINT_1": 0, 
                "POINT_2": 0, 
                "POINT_3": 0, 
                "POINT_4": 0,
                "POINT_5": 0, 
                "POINT_6": 0, 
                "POINT_7": 0, 
                "POINT_8": 0,
                "EDGE_1": 1, 
                "EDGE_2": 1, 
                "EDGE_3": 1, 
                "EDGE_4": 1,
                "EDGE_5": 1, 
                "EDGE_6": 1, 
                "EDGE_7": 1, 
                "EDGE_8": 1,
                "EDGE_9": 1, 
                "EDGE_10": 1, 
                "EDGE_11": 1, 
                "EDGE_12": 1,
                "FACE_1": 2, 
                "FACE_2": 2, 
                "FACE_3": 2, 
                "FACE_4": 2,
                "FACE_5": 2, 
                "FACE_6": 2}

DIM_TO_GROUP = {0: ("POINT_1", "POINT_2", "POINT_3", "POINT_4", "POINT_5",
                    "POINT_6", "POINT_7", "POINT_8"),
                1: ("EDGE_1", "EDGE_2", "EDGE_3", "EDGE_4", "EDGE_5", "EDGE_6",
                    "EDGE_7", "EDGE_8", "EDGE_9", "EDGE_10", "EDGE_11",
                    "EDGE_12"),
                2: ("FACE_1", "FACE_2", "FACE_3", "FACE_4", "FACE_5", 
                    "FACE_6")}

BOUNDARY_MAP = {"POINT_1": [], 
                "POINT_2": [], 
                "POINT_3": [], 
                "POINT_4": [],
                "POINT_5": [], 
                "POINT_6": [], 
                "POINT_7": [], 
                "POINT_8": [],
                "EDGE_1": ["POINT_1", "POINT_2"], 
                "EDGE_2": ["POINT_2", "POINT_3"],
                "EDGE_3": ["POINT_3", "POINT_4"],
                "EDGE_4": ["POINT_4", "POINT_1"],
                "EDGE_5": ["POINT_5", "POINT_6"],
                "EDGE_6": ["POINT_6", "POINT_7"],
                "EDGE_7": ["POINT_7", "POINT_8"],
                "EDGE_9": ["POINT_2", "POINT_6"],
                "EDGE_10": ["POINT_1", "POINT_5"],
                "EDGE_11": ["POINT_3", "POINT_7"],
                "EDGE_12": ["POINT_4", "POINT_8"],
                "FACE_1": ["EDGE_1", "EDGE_2", "EDGE_3", "EDGE_4"],
                "FACE_2": ["EDGE_5", "EDGE_6", "EDGE_7", "EDGE_8"],
                "FACE_3": ["EDGE_10", "EDGE_5", "EDGE_9", "EDGE_1"],
                "FACE_4": ["EDGE_12", "EDGE_7", "EDGE_11", "EDGE_3"],
                "FACE_5": ["EDGE_9", "EDGE_6", "EDGE_11", "EDGE_2"],
                "FACE_6": ["EDGE_10", "EDGE_8", "EDGE_12", "EDGE_4"]}

EL_TYPE_MAP = {}

def initialize_model():
    """
    Initialize the Gmsh API.

    Returns
    -------
    None.

    """
    gmsh.initialize()
    
def finalize_model():
    """
    Finalize the Gmsh API.

    Returns
    -------
    None.

    """
    gmsh.finalize()
    
def show_gmsh_gui(model_name=None):
    """
    Open the Gmsh user interface.

    Parameters
    ----------
    model_name : str, optional
        Name of the model to show in the GUI. The default is None.

    Returns
    -------
    None.

    """
    current = gmsh.model.getCurrent()
    if model_name is not None:
        gmsh.model.setCurrent(model_name)
    gmsh.fltk.run()
    gmsh.model.setCurrent(current)
        
def create_beam_structure(model_name):
    """
    Convert a continuum mesh to a beam mesh.

    Parameters
    ----------
    model_name : str
        Name of the continuum model to convert.

    Returns
    -------
    node_tags : numpy.ndarray
        Node tags of the beam model.
    node_coords : numpy.ndarray
        Coordinates of the node tags of the beam model.
    el_tags : numpy.ndarray
        Element tags of the beam model.
    conn : numpy.ndarray
        Connectivity list of the beam model.

    """
    # TODO: identify element type (e.g. on length of element connectivity)
    
    current = gmsh.model.getCurrent()
    gmsh.model.setCurrent(model_name)
    
    # extract element edges
    gmsh.model.mesh.createEdges()
    edge_tags, edge_nodes = gmsh.model.mesh.getAllEdges()
    edge_nodes = edge_nodes.reshape(-1, 2)
    
    # consider edges as beam elements
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    node_coords = node_coords.reshape(-1, 3)
    el_tags = edge_tags
    conn = edge_nodes
    
    gmsh.model.setCurrent(current)
    return node_tags, node_coords, el_tags, conn

def add_boundary_groups(model_name, group_map):
    """
    Add physical groups for boundary entities to the Gmsh model.

    Parameters
    ----------
    model_name : str
        Name of the Gmsh model.
    group_map : dict
        Map from group names to group entities considering mesh periodicity.

    Returns
    -------
    None.

    """
    current = gmsh.model.getCurrent()
    gmsh.model.setCurrent(model_name)
    
    # make sure all keys are in upper case
    for key in list(group_map.keys()):
        group_map[key.upper()] = group_map.pop(key)
    
    # add physical groups to model
    for name, entities in group_map.items():
        dim = GROUP_TO_DIM[name]
        tag = gmsh.model.addPhysicalGroup(dim, entities)
        gmsh.model.setPhysicalName(dim, tag, name)
    
    gmsh.model.setCurrent(current)

def add_periodic_constraints(model_name, group_map, dx, dy, dz=None):
    """
    Add periodic mesh constraints to the Gmsh model.

    Parameters
    ----------
    model_name : str
        Name of the Gmsh model.
    group_map : dict
        Map from group names to group entities considering mesh periodicity.
    dx : float
        Size of the RVE in x-direction.
    dy : float
        Size of the RVE in y-direction.
    dz : float, optional
        Size of the RVE in z-direction. The default is None.

    Returns
    -------
    None.

    """
    current = gmsh.model.getCurrent()
    gmsh.model.setCurrent(model_name)
    
    # apply periodic constraints on model boundary
    dim = gmsh.model.getDimension()
    if dim == 2:
        transformation = [1, 0, 0, dx, 
                          0, 1, 0, 0, 
                          0, 0, 1, 0, 
                          0, 0, 0, 1]
        gmsh.model.mesh.setPeriodic(1, group_map["EDGE_2"], 
                                    group_map["EDGE_4"], transformation)
        transformation = [1, 0, 0, 0, 
                          0, 1, 0, dy, 
                          0, 0, 1, 0, 
                          0, 0, 0, 1]
        gmsh.model.mesh.setPeriodic(1, group_map["EDGE_3"], 
                                    group_map["EDGE_1"], transformation)
        
    elif dim == 3:
        transformation = [1, 0, 0, dx, 
                          0, 1, 0, 0, 
                          0, 0, 1, 0, 
                          0, 0, 0, 1]
        gmsh.model.mesh.setPeriodic(2, group_map["FACE_2"], 
                                    group_map["FACE_1"], transformation)
        transformation = [1, 0, 0, 0, 
                          0, 1, 0, dy, 
                          0, 0, 1, 0, 
                          0, 0, 0, 1]
        gmsh.model.mesh.setPeriodic(2, group_map["FACE_4"], 
                                    group_map["FACE_3"], transformation)
        transformation = [1, 0, 0, 0, 
                          0, 1, 0, 0, 
                          0, 0, 1, dz, 
                          0, 0, 0, 1]
        gmsh.model.mesh.setPeriodic(2, group_map["FACE_6"], 
                                    group_map["FACE_5"], transformation)
     
    # restore previous model
    gmsh.model.setCurrent(current)
    
def _get_periodic_nodes(dim, tag, include_higher_order_nodes=False):
    """
    Get periodic nodes from physical group information.

    Parameters
    ----------
    dim : int
        Physical group dimension.
    tag : int
        Physical group tag.
    include_higher_order_nodes : bool, optional
        If true, high-order nodes will be in the returned data sets. The 
        default is False.

    Returns
    -------
    tag_master : int
        Tag of the physical master group.
    node_tags : numpy.ndarray
        Node tags in the slave group.
    node_tags_master : numpy.ndarray
        Node_tags in the master group.
    transformation : numpy.ndarray
        Affine transformation matrix of the periodic contraint.

    """
    entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
    
    node_tags = []
    node_tags_master = []
    for e in entities:
        # read periodic contraints
        (tag_master, 
         slave_nodes, 
         master_nodes, 
         transformation) = (gmsh.model.mesh
                            .getPeriodicNodes(dim, e, 
                                              include_higher_order_nodes))
        tag_master = gmsh.model.getPhysicalGroupsForEntity(dim, tag_master)[0]
        node_tags.extend(slave_nodes)
        node_tags_master.extend(master_nodes)
        
    # make sure only unique nodes are read
    _, unique_slave = np.unique(node_tags, return_index=True)
    _, unique_master = np.unique(node_tags_master, return_index=True)
    node_tags = np.array(node_tags)[np.sort(unique_slave)]
    node_tags_master = np.array(node_tags_master)[np.sort(unique_master)]
    
    return tag_master, node_tags, node_tags_master, transformation

def _get_boundary(dim, tag):
    """
    Get boundary groups of a physical group.

    Parameters
    ----------
    dim : int
        Physical group dimension.
    tag : int
        Physical group tag.

    Returns
    -------
    boundary_groups : list
        List of dimension and tags of the boundary groups.

    """
    entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
    entities = [(dim, e) for e in entities]
    boundary = gmsh.model.getBoundary(entities, oriented=False)
    boundary_groups = []
    for e in boundary:
        boundary_group = gmsh.model.getPhysicalGroupsForEntity(*e)
        if len(boundary_group) > 0:
            boundary_groups.append((e[0], boundary_group[0]))
            
    return boundary_groups
        
def create_boundary_sets(model_name, group_map):
    """
    Create node for the application of periodic boundary conditions.

    Parameters
    ----------
    model_name : str
        Name of the Gmsh model.
    group_map : dict
        Map from group names to group entities considering mesh periodicity.

    Returns
    -------
    node_sets : dict
        Map from set names to set nodes considering mesh periodicty.

    """
    def recursively_build_node_sets(group, node_sets, group_map):
        """Recursive generation of node sets."""
        if group[0] == 0:
            return
        
        # get parent group information
        parent_name = gmsh.model.getPhysicalName(*group)
        parent_nodes = node_sets[parent_name]
        
        boundary = _get_boundary(*group)
        to_delete = []
        for g in boundary:
            name = gmsh.model.getPhysicalName(*g)
            
            # collect nodes to be deleted from parent set
            for e in group_map[name]:
                nodes = gmsh.model.mesh.getNodes(g[0], e, includeBoundary=True)
                indices = np.flatnonzero(np.isin(parent_nodes, nodes[0]))
                to_delete.extend(indices)
                
            # add nodes to boundary set
            if name not in node_sets.keys():
                for e in group_map[name]:
                    nodes = gmsh.model.mesh.getNodes(g[0], e, 
                                                     includeBoundary=True)
                    indices = np.flatnonzero(np.isin(parent_nodes, nodes[0]))
                    vals = [node_sets[parent_name][i] for i in indices]
                    node_sets[name] = np.hstack((node_sets.get(name, []), 
                                                  vals))
                
            # call recursive function for child group
            recursively_build_node_sets(g, node_sets, group_map)
        
        # delete boundary nodes from parent set
        node_sets[parent_name] = np.delete(node_sets[parent_name], to_delete)
                
    # activate the requested model
    current = gmsh.model.getCurrent()
    gmsh.model.setCurrent(model_name)
    
    node_sets = {}
    dim = gmsh.model.getDimension()
    groups = gmsh.model.getPhysicalGroups(dim - 1)
    for g in groups:
        # read periodic contraints
        tag_master, slave, master, _ = _get_periodic_nodes(*g)
        if tag_master == g[1]:
            continue
        
        # sort node tags and add them to dictionary
        master_name = gmsh.model.getPhysicalName(dim - 1, tag_master)
        slave_name = gmsh.model.getPhysicalName(*g)
        sorter = np.argsort(master)
        node_sets[master_name] = np.hstack((node_sets.get(master_name, []), 
                                            master[sorter]))
        node_sets[slave_name] = np.hstack((node_sets.get(slave_name, []),
                                           slave[sorter]))

        # proceed with boundary groups
        recursively_build_node_sets(g, node_sets, group_map)
        recursively_build_node_sets((g[0], tag_master), node_sets, group_map)
        
    # make integer arrays
    for name in node_sets.keys():
        node_sets[name] = node_sets[name].astype(int)
    
    gmsh.model.setCurrent(current)
    return node_sets

def write_abq_input(model_name, file_name, abq_el_type, node_sets=None, 
                    element_set=None, ignore_node_sets=None, lattice=False):
    """
    

    Parameters
    ----------
    model_name : TYPE
        DESCRIPTION.
    file_name : TYPE
        DESCRIPTION.
    node_sets : TYPE, optional
        DESCRIPTION. The default is None.
    element_set : TYPE, optional
        DESCRIPTION. The default is None.
    ignore_node_sets : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    path : TYPE
        DESCRIPTION.

    """
    def add_set(file, set_type, set_name, set_items):
        """Add set definitions to the input file."""
        set_type = set_type.upper()
        assert set_type in ("ELSET", "NSET")
        file.write("*{}, {}={}\n".format(set_type, set_type, set_name))
        for i, tag in enumerate(set_items):
            file.write("{}, ".format(tag))
            if (i + 1)%16 == 0:
                file.write("\n")
        file.write("\n")
        
    # get model dimension
    dim = gmsh.model.getDimension()
    
    # TODO: improve lattice handling
    if lattice:
        (node_tags, node_coords, 
         el_tags, conn) = create_beam_structure(model_name)
    
    else:
        # get node information
        node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
        node_coords = node_coords.reshape(-1, 3)
        
        # TODO: make this work also for multiple element types
        # TODO: check if connectivity needs to be reordered for some element types
        # for compatibility between gmsh elements and abaqus elements
        
        # get element information
        el_types, el_tags, conn = gmsh.model.mesh.getElements(dim, 1)
        assert len(el_types) == 1
        el_type = el_types[0]
        el_tags = el_tags[0]
        conn = conn[0]
        
        _, _, _, nen, _, _ = gmsh.model.mesh.getElementProperties(el_type)
        el_tags = el_tags - np.min(el_tags) + 1  # normalize element tags
        conn = conn.reshape(-1, nen)
    
    el_type = abq_el_type  # TODO
    
    # initialize function that checks if elements should be skipped
    ignore_nodes = []
    if ignore_node_sets is not None and node_sets is not None:
        ignore_nodes = np.hstack([node_sets[s] for s in ignore_node_sets])
    ignore = lambda n: np.all(np.in1d(n, ignore_nodes))
    
    # initialize abaqus input file
    file_name, extension = os.path.splitext(file_name)
    file_name += ".inp"
    path = os.path.abspath(file_name)
    with open(path, "w") as file:
        file.write("*PART, NAME=RVE\n")
    
        # add nodes
        file.write("*NODE\n")
        for tag, coords in zip(node_tags, node_coords):
            x, y, z = coords
            file.write("{}, {}, {}, {}\n".format(tag, x, y, z))
            
        # add elements
        file.write("*ELEMENT, TYPE={}\n".format(el_type))
        for tag, nodes in zip(el_tags, conn):
            if ignore(nodes):
                continue
            file.write("{}, ".format(tag))
            for i, n in enumerate(nodes):
                file.write("{}".format(n))
                if i == len(nodes) - 1:
                    file.write("\n")
                    continue
                file.write(", ")
        add_set(file, "NSET", "ALL_NODES", node_tags)
        add_set(file, "ELSET", "ALL_ELEMENTS", el_tags)
        file.write("*END PART\n")
        
    return path









# def _build_node_sets(model_name):
#     """
#     Recursively build the node sets of faces, edges and points.
    
#     The periodic node pairs are conserved by this method.

#     Parameters
#     ----------
#     model_name : str
#         Name of the Gmsh model.

#     Returns
#     -------
#     node_sets : dict
#         Dictionary containing the node sets. Keys are the set names and values
#         are the node tags in the set.

#     """
#     def recursively_build_node_sets(entity, 
#                                     node_sets):
#         """Recursive generation of node sets."""
#         if entity[0] == 0:
#             return
        
#         # get parent entity information
#         parent_name = gmsh.model.getEntityName(*entity)
#         parent_nodes = node_sets[parent_name]
        
#         boundary = gmsh.model.getBoundary([entity], oriented=False)
#         to_delete = []
#         for e in boundary:
#             # get boundary nodes indices in parent
#             nodes = gmsh.model.mesh.getNodes(*e, includeBoundary=True)[0]
#             indices = np.flatnonzero(np.isin(parent_nodes, nodes))
#             to_delete.extend(indices)
            
#             # create entry in dictionary for boundary entity
#             name = gmsh.model.getEntityName(*e)
#             vals = [node_sets[parent_name][i] for i in indices]
#             node_sets[name] = vals
            
#             # call function on child entity
#             recursively_build_node_sets(e, node_sets)
           
#         # remove boundary nodes from parent node set
#         node_sets[parent_name] = np.delete(node_sets[parent_name], to_delete)
        
#     current = gmsh.model.getCurrent()
#     gmsh.model.setCurrent(model_name)
    
#     node_sets = {}
    
#     dim = gmsh.model.getDimension()
#     entities = gmsh.model.getEntities(dim - 1)
#     for e in entities:
#         # read periodic contraints
#         tag_master, slave, master, _ = gmsh.model.mesh.getPeriodicNodes(*e)
#         if tag_master == e[1]:
#             continue
        
#         # sort node tags and add them to dictionary
#         slave_name = gmsh.model.getEntityName(*e)
#         master_name = gmsh.model.getEntityName(e[0], tag_master)
#         sorter = np.argsort(master)
#         node_sets[master_name] = master[sorter].tolist()
#         node_sets[slave_name] = slave[sorter].tolist()
        
#         # print(node_sets)
#         recursively_build_node_sets(e, node_sets)
#         recursively_build_node_sets((e[0], tag_master), node_sets)
    
#     gmsh.model.setCurrent(current)   
#     return node_sets


# def build_lattice_RVE(model_name, dim, edge_length, lc=None, gui=False):
#     """
#     Create of statistical lattice RVE.

#     Parameters
#     ----------
#     model_name : str
#         Name of the model.
#     dim : int
#         Dimension of the model.
#     edge_length : float
#         Edge length of the RVE.
#     lc : float
#         Characteristig cell size.
#     gui : bool, optional
#         If true, the Gmsh GUI will be opened after mesh generation. The 
#         default is False.

#     Returns
#     -------
#     node_tags : numpy.ndarray
#         Tags of the RVE mesh nodes.
#     node_coords : numpy.ndarray
#         Coordinates of the RVE mesh nodes.
#     el_tags : numpy.ndarray
#         Tags of the RVE mesh elements.
#     conn : numpy.ndarray
#         Connectivity list of the RVE mesh.
#     node_sets : dict
#         Node sets of the boundary nodes. The nodes are ordered so that
#         periodic pairs match each other.

#     """
#     # build RVE made out of continuum elements in gmsh
#     if dim == 2:
#         _make_full_RVE_2d(model_name=model_name, gui=gui, 
#                           edge_length=edge_length, lc=lc)
#     else:
#         _make_full_RVE_3d(model_name=model_name, gui=gui, 
#                           edge_length=edge_length, lc=lc)
        
#     # convert continuum mesh to beam mesh
#     node_tags, node_coords, el_tags, conn = _make_beam_structure(model_name)

#     # set up node sets for boundary conditions
#     # node_sets = _build_node_sets(model_name)
#     node_sets = None
    
#     finalize_model()
    
#     return node_tags, node_coords, el_tags, conn, node_sets

