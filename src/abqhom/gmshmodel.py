"""
Module containing the GMshRVE class that is helpful to create a representative
volume element using Gmsh and to prepare the model for homogenization.

"""

import os
import numpy as np
import gmsh
from abqhom.constants import (POINT_GROUPS, EDGE_GROUPS, FACE_GROUPS,
                              PERIODIC_PAIRS, BOUNDARY_CONNECTIVITY,
                              ABAQUS_ELEMENT_TYPES, ABAQUS_CONNECTIVITY_ORDER)


class GmshRVE():
    """
    An RVE representation to prepare a Gmsh model for the Abaqus export.

    Parameters
    ----------
    model_name : str
        The name of the Gmsh model.
    boundary_map : dict[str, tuple[int, list[int]]]
        A map defining the entities being part of the RVE boundary and their
        respective geometric dimensionality

    Attributes
    ----------
    model_name : str
        The name of the Gmsh model.
    node_sets : dict[str, numpy.ndarray[int]]
        A map from set names to set nodes for user-defined physical groups and
        further default node sets.
    element_sets : dict[str, numpy.ndarray[int]]
        A map from set names to set elements for user-defined physical groups
        and further default element sets.
    model_dimension : int
        The geometric dimensionality of the model.
    mesh_dimension : int
        The dimensionality of the mesh of the model

    Example
    -------
    >>> import gmsh
    >>>
    >>> # Build a square-shaped strut-based lattice unit cell in Gmsh
    >>> gmsh.initialize()
    >>> gmsh.model.add("SQUARE_LATTICE")
    >>> p1 = gmsh.model.geo.add_point(0.0, 0.0, 0.0)
    >>> p2 = gmsh.model.geo.add_point(1.0, 0.0, 0.0)
    >>> p3 = gmsh.model.geo.add_point(1.0, 1.0, 0.0)
    >>> p4 = gmsh.model.geo.add_point(0.0, 1.0, 0.0)
    >>> l1 = gmsh.model.geo.add_line(p1, p2)
    >>> l2 = gmsh.model.geo.add_line(p2, p3)
    >>> l3 = gmsh.model.geo.add_line(p3, p4)
    >>> l4 = gmsh.model.geo.add_line(p4, p1)
    >>> gmsh.model.geo.synchronize()
    >>> all_struts_group = gmsh.model.add_physical_group(1, [l1, l2, l3, l4])
    >>> gmsh.model.set_physical_name(1, all_struts_group, "ALL_STRUTS")
    >>>
    >>> # Create a map defining the dimensions and the idetifier tags for the
    >>> # boundary entities
    >>> boundary_map = {"POINT_1": (0, [p1]),
    >>>                 "POINT_2": (0, [p2]),
    >>>                 "POINT_3": (0, [p3]),
    >>>                 "POINT_4": (0, [p4]),
    >>>                 "EDGE_1" : (1, [l1]),
    >>>                 "EDGE_2" : (1, [l2]),
    >>>                 "EDGE_3" : (1, [l3]),
    >>>                 "EDGE_4" : (1, [l4])}
    >>>
    >>> # Create a GmshRVE object and play around with it
    >>> unit_cell = GmshRVE("SQUARE_LATTICE", boudary_map)
    >>> unit_cell.generate_mesh()
    >>> unit_cell.build_node_sets()
    >>> unit_cell.build_element_sets()
    >>> input_file_path = unit_cell.write_to_abaqus_input()

    """

    def __init__(self, model_name, boundary_map):
        """Constructor method for the GmshRVE class."""
        self._model_name = model_name
        self._node_sets = {}
        self._element_sets = {}

        initial_boundary_map = boundary_map.copy()
        for name, (_, entities) in initial_boundary_map.items():
            if len(entities) == 0:
                boundary_map.pop(name)
        self._boundary_map = boundary_map

        if not gmsh.is_initialized():
            gmsh.initialize()

        if model_name not in gmsh.model.list():
            raise ValueError("The Gmsh model {} does not exist.".format(
                model_name))

    @property
    def model_name(self):
        """The name of the Gmsh model."""
        return self._model_name

    @property
    def node_sets(self):
        """Node sets based on physical groups from Gmsh."""
        return self._node_sets

    @property
    def element_sets(self):
        """Element sets based on physical groups from Gmsh."""
        return self._element_sets

    @property
    def model_dimension(self):
        """Get the geometric dimensionality of the RVE model."""

        current = gmsh.model.get_current()
        gmsh.model.set_current(self.model_name)

        if len(gmsh.model.get_entities()) == 0:
            gmsh.model.set_current(current)
            return -1

        # Compute the geometric dimension of the model
        bounding_box = np.array(gmsh.model.get_bounding_box(-1, -1))
        bounding_box = bounding_box.reshape(2, -1)
        delta_x, delta_y, delta_z = bounding_box[1] - bounding_box[0]

        max_delta = max(delta_x, delta_y, delta_z)
        min_deltas = [delta_x, delta_y, delta_z]
        min_deltas.remove(max_delta)
        dim = 3
        for i in min_deltas:
            if i/max_delta < 1e-5:
                dim -= 1

        gmsh.model.set_current(current)
        return dim

    @property
    def mesh_dimension(self):
        """Get the dimensionality of the mesh of the RVE model."""

        current = gmsh.model.get_current()
        gmsh.model.set_current(self.model_name)
        el_types = gmsh.model.mesh.get_elements()[0]
        if len(el_types) == 0:
            gmsh.model.set_current(current)
            return -1
        dim = max((gmsh.model.mesh.get_element_properties(i)[1]
                   for i in el_types))
        gmsh.model.set_current(current)
        return dim

    def __del__(self):
        """Remove the model from Gmsh as soon as the object gets deleted."""

        if gmsh.is_initialized():
            current = gmsh.model.get_current()
            try:
                gmsh.model.set_current(self.model_name)
                gmsh.model.remove()
            except Exception:
                pass

            if current != self.model_name:
                gmsh.model.set_current(current)

    def __str__(self):
        descr = "GmshRVE '{}'\n".format(self.model_name)
        descr += "\tmodel dimension : {}\n".format(self.model_dimension)
        descr += "\tmesh dimension  : {}\n".format(self.mesh_dimension)
        descr += "\tnode sets       : {}\n".format(len(self.node_sets))
        descr += "\telement sets    : {}\n".format(len(self.element_sets))
        return descr

    def activate(self):
        """
        Set the RVE model as the current model in Gmsh.

        Returns
        -------
        None.

        """

        gmsh.model.set_current(self.model_name)

    def show(self):
        """
        Open the RVE model in the Gmsh graphical user interface.

        Returns
        -------
        None.

        """

        current = gmsh.model.get_current()
        gmsh.model.set_current(self.model_name)
        gmsh.fltk.run()
        gmsh.model.set_current(current)

    def build_boundary_map(self):
        raise NotImplementedError("The method is not available yet.")

        dim = self.model_dimension
        entities = gmsh.model.get_entities(dim - 1)
        for e in entities:
            pass

    def generate_mesh(self, element_order=1, periodic=True, algorithm_2d=6,
                      algorithm_3d=1, element_size=None):
        """
        Run the Gmsh meshing routine on the RVE model.

        Parameters
        ----------
        element_order : int
            The order of the generated elements. Possible values are 1 and 2.
            The default is 1.
        periodic : bool, optional
            If True, the mesh will be suitable for application of periodic
            boundary conditions.
        algorithm_2d : int, optional
            The 2d meshing algorithm to use. The available algorithms are
            listed in the Gmsh documentation. The default is 6.
        algorithm_3d : int, optional
            The 3d meshing algorithm to use. The available algorithms are
            listed in the Gmsh documentation. The default is 1.
        element_size : float, optional
            The characteristic element length. If this is None Gmsh will
            set the element size automatically. The default is None.

        Returns
        -------
        None.

        """

        self._add_boundary_groups()
        if periodic:
            self._add_periodic_mesh_constraints()

        # Set element order
        if element_order not in {1, 2}:
            raise ValueError("only first-order and second-order elements are ",
                             "supported")
        gmsh.model.mesh.set_order(element_order)
        gmsh.option.set_number("Mesh.SecondOrderIncomplete", 1)

        # Set meshing algorithms
        gmsh.option.set_number("Mesh.Algorithm", algorithm_2d)
        gmsh.option.set_number("Mesh.Algorithm3D", algorithm_3d)

        # Set mesh size
        if element_size is not None:
            gmsh.option.set_number("Mesh.MeshSizeFromPoints", 0)
            gmsh.option.set_number("Mesh.MeshSizeMin", element_size)
            gmsh.option.set_number("Mesh.MeshSizeMax", element_size)

        mesh_dim = gmsh.model.get_dimension()
        gmsh.model.mesh.generate(mesh_dim)

    def _add_boundary_groups(self):
        """Add physical groups for boundary entities to the Gmsh model."""

        current = gmsh.model.get_current()
        gmsh.model.set_current(self.model_name)

        # Make sure all keys are in upper case
        for key in list(self._boundary_map):
            self._boundary_map[key.upper()] = self._boundary_map.pop(key)

        # Add physical groups to model
        for name, (dim, entities) in self._boundary_map.items():
            tag = gmsh.model.add_physical_group(dim, entities)
            gmsh.model.set_physical_name(dim, tag, name)

        gmsh.model.set_current(current)

    def _add_periodic_mesh_constraints(self):
        """Add periodic mesh constraints to the Gmsh model."""

        def translation(axis, distance):
            """Define the affine matrix for a translation along an axis."""
            transformation = [1, 0, 0, 0,
                              0, 1, 0, 0,
                              0, 0, 1, 0,
                              0, 0, 0, 1]
            transformation[(axis - 1)*4 + 3] = distance
            return transformation

        current = gmsh.model.get_current()
        gmsh.model.set_current(self.model_name)

        # Define transformations
        bounding_box = np.array(
            gmsh.model.get_bounding_box(-1, -1)).reshape(2, -1)
        delta_x, delta_y, delta_z = bounding_box[1] - bounding_box[0]

        transformation_map = {"EDGE_1":  (2, delta_y),
                              "EDGE_3":  (3, delta_z),
                              "EDGE_7":  (2, -delta_y),
                              "EDGE_2":  (1, -delta_x),
                              "EDGE_4":  (3, delta_z),
                              "EDGE_8":  (1, delta_x),
                              "EDGE_9":  (1, delta_x),
                              "EDGE_10": (2, delta_y),
                              "EDGE_11": (1, -delta_x),
                              "FACE_1":  (3, delta_z),
                              "FACE_3":  (2, delta_y),
                              "FACE_5":  (1, delta_x)}

        # Define periodic mesh contraints
        for master, slave in PERIODIC_PAIRS.items():
            if (master not in self._boundary_map
                    or slave not in self._boundary_map):
                continue

            dim = self._boundary_map[master][0]
            if dim not in {1, 2}:
                raise RuntimeError("Periodic mesh constraints are supported "
                                   "for 1d and 2d entities only.")

            gmsh.model.mesh.set_periodic(
                dim, self._boundary_map[slave][1],
                self._boundary_map[master][1],
                translation(*transformation_map[master]))

        gmsh.model.set_current(current)

    def build_node_sets(self):
        """
        Build node sets based on the physical groups from Gmsh.

        The node sets are needed for the correct application of boundary
        conditions.

        Returns
        -------
        None.

        """

        def recursively_assign_nodes(boundary_groups, node_sets):
            """Assign the nodes to the respective boundary groups."""

            for group in boundary_groups:
                dim_tags = self._boundary_map.get(group)
                if dim_tags is None:
                    boundary = BOUNDARY_CONNECTIVITY[group]
                    node_sets = recursively_assign_nodes(boundary, node_sets)
                    continue

                # In case of a point group there are no periodic contraints
                dim, entities = dim_tags
                if group in POINT_GROUPS:
                    entity = entities[0]
                    node = gmsh.model.mesh.get_nodes(dim, entity, True)[0][0]
                    node_sets[(dim, entity)][node] = 0
                    return node_sets

                # In case of lines and faces evaluate the periodic contraints
                # such that the correct order of the nodes is preserved
                for slave_entity in entities:
                    (master_entity,
                     slave_nodes,
                     master_nodes,
                     _) = gmsh.model.mesh.get_periodic_nodes(dim,
                                                             slave_entity,
                                                             True)

                    # If the mesh was created without periodic constraints
                    # the node order is not important anymore
                    if len(master_nodes) == len(slave_nodes) == 0:
                        master_nodes = gmsh.model.mesh.get_nodes(
                            dim, master_entity, True)[0]
                        slave_nodes = gmsh.model.mesh.get_nodes(
                            dim, slave_entity, True)[0]

                    # Assign nodes and their indices in the periodic sets to
                    # the dictionary for storage
                    master_node_set = node_sets[(dim, master_entity)]
                    slave_node_set = node_sets[(dim, slave_entity)]

                    # The master set was already created such that the slave
                    # set will assume the same node order
                    if len(master_node_set) > 0:
                        for master_node, slave_node in zip(master_nodes,
                                                           slave_nodes):
                            slave_node_set[slave_node] = master_node_set[
                                master_node]

                    # The slave set was already created such that the master
                    # set will assume the same node order
                    elif len(slave_node_set) > 0:
                        for master_node, slave_node in zip(master_nodes,
                                                           slave_nodes):
                            master_node_set[master_node] = slave_node_set[
                                slave_node]

                    # Neither the master set nor the slave set exist already
                    else:
                        for index, (master_node, slave_node) in enumerate(
                                zip(master_nodes, slave_nodes)):
                            master_node_set[master_node] = index
                            slave_node_set[slave_node] = index

                # Move on to the boundary of the currently processed group
                boundary = BOUNDARY_CONNECTIVITY[group]
                node_sets = recursively_assign_nodes(boundary, node_sets)

            return node_sets

        def remove_boundary_nodes(parent, children, node_sets):
            """Remove nodes on the boundary of an entity from its node set."""

            for child in children:
                node_sets[parent] = np.setdiff1d(node_sets[parent],
                                                 node_sets[child],
                                                 assume_unique=True)
                node_sets = remove_boundary_nodes(parent,
                                                  BOUNDARY_CONNECTIVITY[child],
                                                  node_sets)
            return node_sets

        def add_physical_groups(node_sets):
            """Add physical groups that were set while creating the model."""

            physical_groups = gmsh.model.get_physical_groups()
            for dim, tag in physical_groups:
                name = gmsh.model.get_physical_name(dim, tag)
                if name not in self._boundary_map:
                    nodes, _ = gmsh.model.mesh.get_nodes_for_physical_group(
                        dim, tag)
                    node_sets[name] = np.array(nodes)
            return node_sets

        # Activate the model
        current = gmsh.model.get_current()
        gmsh.model.set_current(self.model_name)

        # Return if there are no nodes in the mesh
        if self.mesh_dimension < 0:
            return

        # Recursively create the node sets starting with the highest dimension
        node_sets = {(dim, tag): {} for (dim, entities) in
                     self._boundary_map.values() for tag in entities}
        if self.model_dimension == 3:
            node_sets = recursively_assign_nodes(FACE_GROUPS, node_sets)
        elif self.model_dimension == 2:
            node_sets = recursively_assign_nodes(EDGE_GROUPS, node_sets)
        else:
            raise RuntimeError("Only 2d and 3d models are supported.")

        # Store node sets in arrays such that periodic pairs are preserved and
        # make sure nodes appear only once per set
        node_sets_per_group = {}
        for group, (dim, entities) in self._boundary_map.items():
            nodes = []
            for entity in entities:
                sorter = np.argsort(node_sets[(dim, entity)].values())
                nodes.extend(np.array(
                    list(node_sets[(dim, entity)].keys()))[sorter])
            _, indices = np.unique(nodes, return_index=True)
            node_sets_per_group[group] = np.array(nodes)[np.sort(indices)]
        node_sets = node_sets_per_group

        # Remove nodes being part of the boundary of a group
        for group in self._boundary_map:
            node_sets = remove_boundary_nodes(group,
                                              BOUNDARY_CONNECTIVITY[group],
                                              node_sets)

        # Add additional sets collecting nodes for points, edges and faces
        boundary_groups = POINT_GROUPS + EDGE_GROUPS + FACE_GROUPS
        node_sets = self._combine_existing_sets(
            node_sets,
            ["POINT_NODES", "EDGE_NODES", "FACE_NODES", "BOUNDARY_NODES"],
            [POINT_GROUPS, EDGE_GROUPS, FACE_GROUPS, boundary_groups])

        # Mark boundary nodes and inner nodes
        # Loop over entities instead of directly getting all nodes to avoid
        # selecting nodes on orphan points, e.g. control points of splines
        entities = gmsh.model.get_entities(gmsh.model.get_dimension())
        all_nodes = []
        for dim, entity in entities:
            all_nodes.extend(gmsh.model.mesh.get_nodes(dim, entity, True)[0])
        node_sets["ALL_NODES"] = np.array(all_nodes, dtype=int)
        node_sets["INNER_NODES"] = np.setdiff1d(node_sets["ALL_NODES"],
                                                node_sets["BOUNDARY_NODES"])

        # Add further physical groups that were set while creating the model
        node_sets = add_physical_groups(node_sets)

        # Assign sets to the object property and make sure nodes are integers
        for name in node_sets.keys():
            self.node_sets[name] = np.array(node_sets[name], dtype=int)

        # Restore previous model
        gmsh.model.set_current(current)

    def build_element_sets(self):
        """
        Build node sets based on the physical groups from Gmsh.

        The element sets are needed to assign e.g. material to the elements.

        Returns
        -------
        None.

        """

        # Activate the model
        current = gmsh.model.get_current()
        gmsh.model.set_current(self.model_name)

        # Return if there are no elements in the mesh
        if self.mesh_dimension < 0:
            return

        # Find all groups that need to be processed based on the mesh dimension
        relevant_groups = gmsh.model.get_physical_groups(self.mesh_dimension)

        # Initialize element sets
        element_sets = {gmsh.model.get_physical_name(*group): []
                        for group in relevant_groups}
        for group in self._boundary_map:
            element_sets[group] = []

        # Add all sets containing elements of the highest available dimension
        for group in relevant_groups:
            name = gmsh.model.get_physical_name(*group)
            entities = gmsh.model.get_entities_for_physical_group(*group)
            for entity in entities:
                elements = gmsh.model.mesh.get_elements(group[0], entity)[1]
                elements = [e for el_type in elements for e in el_type]
                element_sets[name].extend(elements)
            element_sets[name] = np.unique(element_sets[name])

        # Add additional sets collecting elements for points, edges and faces
        boundary_groups = POINT_GROUPS + EDGE_GROUPS + FACE_GROUPS
        element_sets = self._combine_existing_sets(
            element_sets,
            ["POINT_ELEMENTS", "EDGE_ELEMENTS", "FACE_ELEMENTS",
             "BOUNDARY_ELEMENTS"],
            [POINT_GROUPS, EDGE_GROUPS, FACE_GROUPS, boundary_groups])

        # Mark boundary nodes and inner nodes
        all_elements = gmsh.model.mesh.get_elements(self.mesh_dimension)[1]
        all_elements = np.array([e for el_type in all_elements
                                 for e in el_type])
        element_sets["ALL_ELEMENTS"] = all_elements
        element_sets["INNER_ELEMENTS"] = np.setdiff1d(
            all_elements, element_sets["BOUNDARY_ELEMENTS"])

        # Assign sets to the object property and make sure elements are
        # integers
        for name in element_sets.keys():
            self.element_sets[name] = np.array(element_sets[name], dtype=int)

        # Restore previous model
        gmsh.model.setCurrent(current)

    def _combine_existing_sets(self, set_dictionary, set_names,
                               set_collections):
        """Add node or element sets combining multiple existing sets."""

        for set_collection, set_name in zip(
                set_collections, set_names):
            set_items = []
            for i in set_collection:
                items = set_dictionary.get(i)
                if items is not None:
                    set_items.extend(items)

            set_dictionary[set_name] = np.array(set_items, dtype=int)

        return set_dictionary

    def write_to_abaqus_input(self, el_type_map=None):
        """
        Export the Gmsh mesh to Abaqus.

        Node and element sets will not be written to the input file and need to
        be added to the model subsequently directly within Abaqus.

        Parameters
        ----------
        el_type_map : dict[int, str], optional
            Map from Gmsh element types to Abaqus element types. If this is
            None a default map will be used. The default is None.

        Returns
        -------
        path : str
            Absolute path to the generated Abaqus input file.

        """

        # Activate the model
        current = gmsh.model.get_current()
        gmsh.model.set_current(self.model_name)

        # Load default element type map
        if el_type_map is None:
            el_type_map = ABAQUS_ELEMENT_TYPES

        # Get node information
        # Loop over entities instead of directly getting all nodes to avoid
        # selecting nodes on orphan points, e.g. control points of splines
        entities = gmsh.model.get_entities(gmsh.model.get_dimension())
        node_tags, node_coords = [], []
        for dim, entity in entities:
            entity_nodes, entity_node_coords, _ = gmsh.model.mesh.get_nodes(
                dim, entity, True)
            node_tags.extend(entity_nodes)
            node_coords.extend(entity_node_coords)
        node_tags = np.array(node_tags, dtype=int)
        node_coords = np.array(node_coords, dtype=float).reshape(-1, 3)

        # Get element information
        el_types, el_tags, conn = gmsh.model.mesh.getElements(
            self.mesh_dimension)

        # Reorder the connectivity list in terms of the Abaqus node order
        for i, el_type in enumerate(el_types):
            _, _, _, nen, _, _ = gmsh.model.mesh.getElementProperties(el_type)
            conn[i] = np.array(conn).reshape(-1, nen)
            conn[i] = conn[i][:, ABAQUS_CONNECTIVITY_ORDER[el_type]]

        # Initialize Abaqus input file
        file_name = self.model_name + ".inp"
        path = os.path.abspath(file_name)
        with open(path, "w") as file:
            file.write("*PART, NAME={}\n".format(self.model_name))

            # Add nodes
            file.write("*NODE\n")
            for tag, coords in zip(node_tags, node_coords):
                x, y, z = coords
                file.write("{}, {}, {}, {}\n".format(tag, x, y, z))

            # Add elements
            for el_type, type_tags, type_conn in zip(el_types, el_tags, conn):
                file.write("*ELEMENT, TYPE={}\n".format(el_type_map[el_type]))
                for tag, nodes in zip(type_tags, type_conn):
                    file.write("{}".format(tag))
                    for i, node in enumerate(nodes):
                        file.write(", {}".format(node))
                    file.write("\n")

            file.write("*END PART\n")

        # Restore the previous model
        gmsh.model.set_current(current)

        return path

    def write_to_matlab_input(self):
        """
        Export the entities of the Gmsh model to an input file for the
        voxel-based Matlab solver "3D Homogenization of cellular materials"
        written by Guoying Dong.

        The solver is available from the Matlab File Exchange:
        https://de.mathworks.com/matlabcentral/fileexchange/67457-3d-homogenization-of-cellular-materials

        Returns
        -------
        path : str
            Absolute path to the Matlab input file.

        """

        # Activate the model
        current = gmsh.model.get_current()
        gmsh.model.set_current(self.model_name)

        if not self.model_dimension == 3 and self.mesh_dimension == 1:
            raise RuntimeError("Matlab support is only available for 3d "
                               "lattice unit cells.")

        # Compute RVE size
        bounding_box = np.array(
            gmsh.model.get_bounding_box(-1, -1)).reshape(2, -1)
        delta_x, delta_y, delta_z = bounding_box[1] - bounding_box[0]
        x_min = bounding_box[0][0]
        y_min = bounding_box[0][1]
        z_min = bounding_box[0][2]

        # Get grid point information
        grid = gmsh.model.get_entities(0)
        grid_coords = []
        for dim, tag in grid:
            grid_coords.append(gmsh.model.get_value(dim, tag, []))
        grid_tag_map = {tag: i + 1 for i, (_, tag) in enumerate(grid)}
        grid_tags = [i + 1 for i in range(len(grid))]

        # Normalize grid coordinates
        for coords in grid_coords:
            coords[0] = (coords[0] - x_min)/delta_x
            coords[1] = (coords[1] - y_min)/delta_y
            coords[2] = (coords[2] - z_min)/delta_z

        # Get connectivity information
        struts = gmsh.model.get_entities(1)
        connectivity = []
        for dim, tag in struts:
            _, nodes = gmsh.model.get_adjacencies(dim, tag)
            assert len(nodes) == 2
            nodes = [grid_tag_map[nodes[0]], grid_tag_map[nodes[1]]]
            connectivity.append(nodes)
        strut_tags = [i + 1 for i in range(len(struts))]

        # Initialize output file
        file_name = self.model_name + ".txt"
        path = os.path.abspath(file_name)
        with open(path, "w") as file:

            # Add grid nodes
            file.write("//Grid\tID\tx\ty\tz\n")
            for tag, coords in zip(grid_tags, grid_coords):
                x, y, z = coords
                file.write("GRID\t{}\t{}\t{}\t{}\n".format(tag, x, y, z))

            # Add grid edges
            file.write("//Strut\tID\tStart\tEnd\n")
            for tag, conn in zip(strut_tags, connectivity):
                file.write("STRUT\t{}\t{}\t{}\n".format(tag, conn[0], conn[1]))

        # Restore previous model
        gmsh.model.set_current(current)

        return path
