import numpy as np
import gmsh
import os

def _make_full_RVE_2d(model_name="RVE2d", 
                 edge_length=1.0, 
                 lc=None, 
                 gui=False):
    """
    Create a 2d RVE.
    
    The RVE is meshed with continuum first-order triangular elements.

    Parameters
    ----------
    modelname : str, optional
        Name of the Gmsh model. The default is "RVE2d".
    edge_length : float, optional
        Edge length of the square RVE. Default is 1.0.
    lc : float, optional
        Characteristic element length. Default is None.
    gui : bool, optional
        If true, the Gmsh user interface will be started. The default is False.

    Returns
    -------
    modelname : str
        Name of the Gmsh model.

    """
    # initialize gmsh and add a new model
    if not gmsh.isInitialized():
        gmsh.initialize()
    gmsh.model.add(model_name)
    
    # RVE will be a unit box
    gmsh.model.occ.addRectangle(0.0, 0.0, 0.0, edge_length, edge_length)
    gmsh.model.occ.synchronize()
    
    # set entity names
    _set_entity_names(0, "POINT")
    _set_entity_names(1, "EDGE")
    _set_entity_names(2, "FACE")
    
    # add periodic boundary constraints
    transformation = [1, 0, 0, edge_length, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    gmsh.model.mesh.setPeriodic(1, [2], [4], transformation)
    transformation = [1, 0, 0, 0, 0, 1, 0, edge_length, 0, 0, 1, 0, 0, 0, 0, 1]
    gmsh.model.mesh.setPeriodic(1, [3], [1], transformation)
    
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
    
    return model_name, edge_length
    
def _make_full_RVE_3d(model_name="RVE3d", 
                      edge_length=1.0, 
                      lc=None, 
                      gui=False):
    """
    Create a 3d RVE.
    
    The RVE is meshed with continuum first-order tetrahedral elements.

    Parameters
    ----------
    modelname : str, optional
        Name of the Gmsh model. The default is "RVE3d".
    edge_length : float, optional
        Edge length of the cube RVE. Default is 1.0.
    lc : float
        Characteristic element length. Default is None.
    gui : bool, optional
        If true, the Gmsh user interface will be started. The default is False.

    Returns
    -------
    modelname : str
        Name of the Gmsh model.

    """
    # initialize gmsh and add a new model
    if not gmsh.isInitialized():
        gmsh.initialize()
    gmsh.model.add(model_name)

    # RVE will be a unit box
    gmsh.model.occ.addBox(0.0, 0.0, 0.0, edge_length, edge_length, edge_length)
    gmsh.model.occ.synchronize()

    # set entity names
    _set_entity_names(0, "POINT")
    _set_entity_names(1, "EDGE")
    _set_entity_names(2, "FACE")
    _set_entity_names(3, "VOLUME")
    
    # add periodic boundary constraints
    transformation = [1, 0, 0, edge_length, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    gmsh.model.mesh.setPeriodic(2, [2], [1], transformation)
    transformation = [1, 0, 0, 0, 0, 1, 0, edge_length, 0, 0, 1, 0, 0, 0, 0, 1]
    gmsh.model.mesh.setPeriodic(2, [4], [3], transformation)
    transformation = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, edge_length, 0, 0, 0, 1]
    gmsh.model.mesh.setPeriodic(2, [6], [5], transformation)
    
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

    return model_name, edge_length

def _make_hole_RVE_2d(model_name="RVE2d", 
                      edge_length=1.0, 
                      radius=0.25,
                      lc=None, 
                      gui=False):
    """
    Create a 2d RVE with a hole in the middle.
    
    The RVE is meshed with continuum first-order triangular elements.

    Parameters
    ----------
    modelname : str, optional
        Name of the Gmsh model. The default is "RVE2d".
    edge_length : float, optional
        Edge length of the square RVE. Default is 1.0.
    radius : float, optional
        Radius of the hole. Default is 0.25.
    lc : float, optional
        Characteristic element length. Default is None.
    gui : bool, optional
        If true, the Gmsh user interface will be started. The default is False.

    Returns
    -------
    modelname : str
        Name of the Gmsh model.

    """
    # initialize gmsh and add a new model
    if not gmsh.isInitialized():
        gmsh.initialize()
    gmsh.model.add(model_name)
    
    # RVE will be a unit box
    gmsh.model.occ.addRectangle(0.0, 0.0, 0.0, edge_length, edge_length)
    gmsh.model.occ.addDisk(0.5*edge_length, 0.5*edge_length, 0.0,
                           radius, radius)
    gmsh.model.occ.cut([(2, 1)], [(2, 2)])
    gmsh.model.occ.synchronize()
    
    # set entity names
    _set_entity_names(0, "POINT")
    _set_entity_names(1, "EDGE")
    _set_entity_names(2, "FACE")
    
    # add periodic boundary constraints
    transformation = [1, 0, 0, edge_length, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    gmsh.model.mesh.setPeriodic(1, [3], [2], transformation)
    transformation = [1, 0, 0, 0, 0, 1, 0, edge_length, 0, 0, 1, 0, 0, 0, 0, 1]
    gmsh.model.mesh.setPeriodic(1, [4], [1], transformation)
    
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
    
    return model_name, edge_length

def _make_hole_RVE_3d(model_name="RVE3d", 
                      edge_length=1.0, 
                      radius=0.25,
                      lc=None, 
                      gui=False):
    """
    Create a 3d RVE with a hole in the middle.
    
    The RVE is meshed with continuum first-order tetrahedral elements.

    Parameters
    ----------
    modelname : str, optional
        Name of the Gmsh model. The default is "RVE2d".
    edge_length : float, optional
        Edge length of the cube RVE. Default is 1.0.
    radius : float, optional
        Radius of the hole. Default is 0.25.
    lc : float, optional
        Characteristic element length. Default is None.
    gui : bool, optional
        If true, the Gmsh user interface will be started. The default is False.

    Returns
    -------
    modelname : str
        Name of the Gmsh model.

    """
    # initialize gmsh and add a new model
    if not gmsh.isInitialized():
        gmsh.initialize()
    gmsh.model.add(model_name)

    # RVE will be a unit box
    gmsh.model.occ.addBox(0.0, 0.0, 0.0, edge_length, edge_length, edge_length)
    gmsh.model.occ.addSphere(0.5*edge_length, 0.5*edge_length, 0.5*edge_length,
                             radius)
    gmsh.model.occ.cut([(3, 1)], [(3, 2)])
    gmsh.model.occ.synchronize()

    # set entity names
    _set_entity_names(0, "POINT")
    _set_entity_names(1, "EDGE")
    _set_entity_names(2, "FACE")
    _set_entity_names(3, "VOLUME")
    
    # add periodic boundary constraints
    transformation = [1, 0, 0, edge_length, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    gmsh.model.mesh.setPeriodic(2, [6], [1], transformation)
    transformation = [1, 0, 0, 0, 0, 1, 0, edge_length, 0, 0, 1, 0, 0, 0, 0, 1]
    gmsh.model.mesh.setPeriodic(2, [4], [2], transformation)
    transformation = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, edge_length, 0, 0, 0, 1]
    gmsh.model.mesh.setPeriodic(2, [3], [5], transformation)
    
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

    return model_name, edge_length

def _set_entity_names(dim, 
                      base_name=""):
    """
    Set the name for all entities of dimension `dim`.
    
    `base_name` will be used as a prefix.

    Parameters
    ----------
    dim : int
        Dimension of the entities.
    base_name : str
        Prefix of each entitie's name.

    Returns
    -------
    None.

    """
    for i, e in enumerate(gmsh.model.getEntities(dim)):
        gmsh.model.setEntityName(e[0], e[1], base_name + "_{}".format(i + 1))
        
def _make_beam_structure(model_name):
    """
    Convert the continuum mesh to a beam mesh.

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

def _build_node_sets(model_name):
    """
    Recursively build the node sets of faces, edges and points.
    
    The periodic node pairs are conserved by this method.

    Parameters
    ----------
    model_name : str
        Name of the Gmsh model.

    Returns
    -------
    node_sets : dict
        Dictionary containing the node sets. Keys are the set names and values
        are the node tags in the set.

    """
    def recursively_build_node_sets(entity, 
                                    node_sets):
        """Recursive generation of node sets."""
        if entity[0] == 0:
            return
        
        # get parent entity information
        parent_name = gmsh.model.getEntityName(*entity)
        parent_nodes = node_sets[parent_name]
        
        boundary = gmsh.model.getBoundary([entity], oriented=False)
        to_delete = []
        for e in boundary:
            # get boundary nodes indices in parent
            nodes = gmsh.model.mesh.getNodes(*e, includeBoundary=True)[0]
            indices = np.flatnonzero(np.isin(parent_nodes, nodes))
            to_delete.extend(indices)
            
            # create entry in dictionary for boundary entity
            name = gmsh.model.getEntityName(*e)
            vals = [node_sets[parent_name][i] for i in indices]
            node_sets[name] = vals
            
            # call function on child entity
            recursively_build_node_sets(e, node_sets)
           
        # remove boundary nodes from parent node set
        node_sets[parent_name] = np.delete(node_sets[parent_name], to_delete)
        
    current = gmsh.model.getCurrent()
    gmsh.model.setCurrent(model_name)
    
    node_sets = {}
    
    dim = gmsh.model.getDimension()
    entities = gmsh.model.getEntities(dim - 1)
    for e in entities:
        # read periodic contraints
        tag_master, slave, master, _ = gmsh.model.mesh.getPeriodicNodes(*e)
        if tag_master == e[1]:
            continue
        
        # sort node tags and add them to dictionary
        slave_name = gmsh.model.getEntityName(*e)
        master_name = gmsh.model.getEntityName(e[0], tag_master)
        sorter = np.argsort(master)
        node_sets[master_name] = master[sorter].tolist()
        node_sets[slave_name] = slave[sorter].tolist()
        
        # print(node_sets)
        recursively_build_node_sets(e, node_sets)
        recursively_build_node_sets((e[0], tag_master), node_sets)
    
    gmsh.model.setCurrent(current)   
    return node_sets

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

def build_lattice_RVE(model_name, dim, edge_length, lc=None, gui=False):
    """
    Create of statistical lattice RVE.

    Parameters
    ----------
    model_name : str
        Name of the model.
    dim : int
        Dimension of the model.
    edge_length : float
        Edge length of the RVE.
    lc : float
        Characteristig cell size.
    gui : bool, optional
        If true, the Gmsh GUI will be opened after mesh generation. The 
        default is False.

    Returns
    -------
    node_tags : numpy.ndarray
        Tags of the RVE mesh nodes.
    node_coords : numpy.ndarray
        Coordinates of the RVE mesh nodes.
    el_tags : numpy.ndarray
        Tags of the RVE mesh elements.
    conn : numpy.ndarray
        Connectivity list of the RVE mesh.
    node_sets : dict
        Node sets of the boundary nodes. The nodes are ordered so that
        periodic pairs match each other.

    """
    # build RVE made out of continuum elements in gmsh
    if dim == 2:
        _make_full_RVE_2d(model_name=model_name, gui=gui, 
                          edge_length=edge_length, lc=lc)
    else:
        _make_full_RVE_3d(model_name=model_name, gui=gui, 
                          edge_length=edge_length, lc=lc)
        
    # convert continuum mesh to beam mesh
    node_tags, node_coords, el_tags, conn = _make_beam_structure(model_name)

    # set up node sets for boundary conditions
    node_sets = _build_node_sets(model_name)
    
    finalize_model()
    
    return node_tags, node_coords, el_tags, conn, node_sets

def build_hole_RVE(model_name, dim, edge_length, radius, lc=None, gui=False):
    """
    Create a RVE with a hole in its center.

    Parameters
    ----------
    model_name : str
        Name of the model.
    dim : int
        Dimension of the model.
    edge_length : float
        Edge length of the RVE.
    radius : float
        Radius of the hole.
    lc : float
        Characteristig cell size.
    gui : bool, optional
        If true, the Gmsh GUI will be opened after mesh generation. The 
        default is False.

    Returns
    -------
    node_tags : numpy.ndarray
        Tags of the RVE mesh nodes.
    node_coords : numpy.ndarray
        Coordinates of the RVE mesh nodes.
    el_tags : numpy.ndarray
        Tags of the RVE mesh elements.
    conn : numpy.ndarray
        Connectivity list of the RVE mesh.
    node_sets : dict
        Node sets of the boundary nodes. The nodes are ordered so that
        periodic pairs match each other.

    """
    # build RVE made out of continuum elements in gmsh
    if dim == 2:
        _make_hole_RVE_2d(model_name=model_name, gui=gui, 
                          edge_length=edge_length, radius=radius, lc=lc)
    else:
        _make_full_RVE_3d(model_name=model_name, gui=gui, 
                          edge_length=edge_length, radius=radius, lc=lc)
    
    # extract element edges
    gmsh.model.mesh.createEdges()
    _, el_tags, conn = gmsh.model.mesh.getElements(dim, 1)
    el_tags = el_tags[0] - np.min(el_tags[0]) + 1
    conn = conn[0].reshape(-1, 3 if dim == 2 else 4)
    
    # extract nodes
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    node_coords = node_coords.reshape(-1, 3)
    
    # set up node sets for boundary conditions
    node_sets = _build_node_sets(model_name)
    
    finalize_model()
    
    return node_tags, node_coords, el_tags, conn, node_sets

def RVE_to_abaqus_input(size, node_tags, node_coords, el_tags, conn, el_type, 
                        file):
    """
    Export an RVE mesh representation to Abaqus input file.

    Parameters
    ----------
    size : float
        Edge length of the RVE.
    node_tags : numpy.ndarray
        Tags of the RVE mesh nodes.
    node_coords : numpy.ndarray
        Coordinates of the RVE mesh nodes.
    el_tags : numpy.ndarray
        Tags of the RVE mesh elements.
    conn : numpy.ndarray
        Connectivity list of the RVE mesh.
    el_type : str
        Abaqus element identifier, e.g. "CPE3".
    file : str
        Name of the Abaqus input file.

    Returns
    -------
    path : str
        Absolute path of the Abaqus input file.

    """
    def add_set(content, set_type, set_name, set_items):
        """Add set definitions to the input file."""
        if set_type.upper() not in ("NSET", "ELSET"):
            raise ValueError("{} set is not supported.".format(set_type))
            
        content.write(
            "*{}, {}={}\n".format(set_type, set_type, set_name))
        for i, tag in enumerate(set_items):
            content.write("{}, ".format(tag))
            if (i + 1)%16 == 0:
                content.write("\n")
        content.write("\n")
        
    def add_dummy_node(content, name, tag, coords):
        """Add dummy node for periodic BCs to the input file."""
        x, y, z = coords
        content.write("*PART, NAME={}\n".format(name))
        content.write("*NODE\n")
        content.write("{}, {}, {}, {}\n".format(tag, x, y, z))
        add_set(content, "NSET", "BC", [tag])
        content.write("*END PART\n")
        
    # initialize abaqus input file
    file, extension = os.path.splitext(file)
    file += ".inp"
    path = os.path.abspath(file)
    
    max_node_tag = int(np.max(node_tags))
    with open(path, "w") as content:
        content.write("*PART, NAME=RVE\n")
    
        # add nodes
        content.write("*NODE\n")
        for tag, coords in zip(node_tags, node_coords):
            x, y, z = coords
            content.write("{}, {}, {}, {}\n".format(tag, x, y, z))
            
        # add elements
        content.write("*ELEMENT, TYPE={}\n".format(el_type))
        for tag, nodes in zip(el_tags, conn):
            content.write("{}, ".format(tag))
            for i, n in enumerate(nodes):
                content.write("{}".format(n))
                if i == len(nodes) - 1:
                    content.write("\n")
                    continue
                content.write(", ")
        add_set(content, "NSET", "ALL_NODES", node_tags)
        add_set(content, "ELSET", "ALL_ELEMENTS", el_tags)
        content.write("*END PART\n")
        
        # add dummy nodes to apply boundary conditions
        tag = max_node_tag + 1
        coords = [0.5*size, -1.0, 0.0]
        add_dummy_node(content, "DUMMY_1", tag, coords)
        
        tag = max_node_tag + 2
        coords = [-1.0, 0.5*size, 0.0]
        add_dummy_node(content, "DUMMY_2", tag, coords)
        
        tag = max_node_tag + 3
        coords = [-1.0, 0.0, 0.5*size]
        add_dummy_node(content, "DUMMY_3", tag, coords)
        
    return path