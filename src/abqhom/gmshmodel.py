"""
Module containing the GMshRVE class that is helpful to create a representative
volume element using Gmsh and to prepare the model for homogenization.

"""

import os
from collections import UserDict
import numpy as np
import gmsh
from abqhom.constants import (POINT_GROUPS, EDGE_GROUPS, FACE_GROUPS,
                              PERIODIC_PAIRS, ABAQUS_ELEMENT_TYPES,
                              ABAQUS_CONNECTIVITY_ORDER, BOUNDARY_CONNECTIVITY,
                              SLAVE_GROUPS_OFFSET_AXES)


class GmshRVE():
    """
    An RVE representation to prepare a Gmsh model for the Abaqus export.

    The RVE must be rectangular or box-shaped and must have line or surface
    entities connected to its corner vertices.

    Parameters
    ----------
    model_name : str
        The name of the underlying Gmsh model.
    boundary_map : Boundary, optional
        A map defining the entities being part of the RVE boundary groups.
        The key is the group name and the value is a list of Gmsh dimtags
        of the entities being part of the group. If boundary_map is None,
        abqhom will generate a mapping automatically. The default is None.

    Attributes
    ----------
    model_name : str
        The name of the underlying Gmsh model.
    boundary : Boundary
        The boundary definition of the RVE.
    periodic_boundary : bool
        A flag to mark whether the model is periodic.
    periodic_mesh : bool
        A flag to mark whether the model has periodic mesh constraints.
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
    >>> from abqhom.gmshmodel import GmshRVE
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
    >>> # Create a GmshRVE object and play around with it
    >>> unit_cell = GmshRVE("SQUARE_LATTICE")
    >>> unit_cell.add_periodic_mesh_constraints()
    >>> unit_cell.generate_mesh()
    >>> unit_cell.build_node_sets()
    >>> unit_cell.build_element_sets()
    >>> input_file_path = unit_cell.write_abaqus_input()

    """

    def __init__(self, model_name, boundary=None):
        if not gmsh.is_initialized():
            gmsh.initialize()
        if model_name not in gmsh.model.list():
            raise ValueError("The Gmsh model {} does not exist.".format(
                model_name))

        # Assign attributes
        self._node_sets = {}
        self._element_sets = {}
        self._boundary = boundary
        self._periodic_mesh = False

        # Create the boundary map
        if boundary is None:
            self._boundary = Boundary(model_name)

        # Check if the model has entities of at least dimension 1
        self.activate()
        if len(gmsh.model.get_entities(1)) == 0:
            raise RuntimeError(
                "The model must contain entities of at least dimension 1.")

        # Check if the model is at least a 2d model
        if self.model_dimension < 2:
            raise RuntimeError("The model must be of at least dimension 2.")

        # Make sure the point groups are not empty
        points = POINT_GROUPS
        if self.model_dimension == 2:
            points = ["POINT_1", "POINT_2", "POINT_3", "POINT_4"]
        for p in points:
            if len(self.boundary[p]) == 0:
                raise RuntimeError(
                    "Point groups {} must not be empty.".format(set(points)))

        # Add physical boundary groups to the model
        self._add_boundary_groups()

    @property
    def model_name(self):
        """The name of the Gmsh model."""
        return self._boundary.model_name

    @property
    def boundary(self):
        """The boundary definition of the RVE."""
        return self._boundary

    @property
    def periodic_boundary(self):
        """Flag to mark whether the model is periodic."""
        return self._boundary.periodic

    @property
    def periodic_mesh(self):
        """Flag to mark whether the model has periodic mesh constraints."""
        return self._periodic_mesh

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
        self.activate()

        if len(gmsh.model.get_entities()) == 0:
            return -1

        # Compute the geometric dimension of the model
        tol = gmsh.option.get_number("Geometry.Tolerance")
        bounding_box = np.array(self.get_bounding_box_from_points())
        bounding_box = bounding_box.reshape(2, -1)
        delta_x, delta_y, delta_z = bounding_box[1] - bounding_box[0]
        dim = 3 - np.isclose([delta_x, delta_y, delta_z], 0.0, atol=tol).sum()

        return dim

    @property
    def mesh_dimension(self):
        """Get the dimensionality of the mesh of the RVE model."""

        self.activate()
        el_types = gmsh.model.mesh.get_elements()[0]
        if len(el_types) == 0:
            return -1
        dim = max((gmsh.model.mesh.get_element_properties(i)[1]
                   for i in el_types))
        return dim

    def __str__(self):
        descr = "GmshRVE '{}'\n".format(self.model_name)
        descr += "\tmodel dimension : {}\n".format(self.model_dimension)
        descr += "\tmesh dimension  : {}\n".format(self.mesh_dimension)
        descr += "\tnode sets       : {}\n".format(len(self.node_sets))
        descr += "\telement sets    : {}".format(len(self.element_sets))
        return descr

    def __repr__(self):
        return self.__str__()

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

        self.activate()
        gmsh.fltk.run()

    def get_bounding_box_from_points(self):
        """
        Compute the bounding box of the geometric CAD model from points.

        Returns
        -------
        x_min : float
            The minimum x-value of the bounding box.
        y_min : float
            The minimum y-value of the bounding box.
        z_min : float
            The minimum z-value of the bounding box.
        x_max : float
            The maximum x-value of the bounding box.
        y_max : float
            The maximum y-value of the bounding box.
        z_max : float
            The maximum z-value of the bounding box.

        """

        (x_min, y_min, z_min,
         x_max, y_max, z_max) = self.boundary.get_bounding_box_from_points()
        return x_min, y_min, z_min, x_max, y_max, z_max

    def _add_boundary_groups(self):
        """Add physical groups for boundary entities to the Gmsh model."""

        # Add physical groups to model
        self.activate()
        for name, dimtags in self.boundary.items():
            for dim in range(3):
                entities = [j for i, j in dimtags if i == dim]
                tag = gmsh.model.add_physical_group(dim, entities)
                gmsh.model.set_physical_name(dim, tag, name)

    def generate_mesh(self, element_order=1, algorithm_2d=6, algorithm_3d=1,
                      element_size=None):
        """
        Run the Gmsh meshing routine on the RVE model.

        Parameters
        ----------
        element_order : int
            The order of the generated elements. Possible values are 1 and 2.
            The default is 1.
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

        # Set element order
        if element_order not in {1, 2}:
            raise ValueError("Only first-order and second-order elements are "
                             "supported.")
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

    def add_periodic_mesh_constraints(self):
        """
        Add periodic mesh constraints to the Gmsh model.

        Returns
        -------
        None.

        """

        def translation(axes, edge_lengths):
            """Define the affine matrix for a translation along an axis."""
            transformation = [1, 0, 0, 0,
                              0, 1, 0, 0,
                              0, 0, 1, 0,
                              0, 0, 0, 1]
            for a in axes:
                transformation[a*4 + 3] = edge_lengths[a]
            return transformation

        if not self.periodic_boundary:
            raise RuntimeError(
                "The model does not have a periodic boundary formulation.")

        # Define transformations
        bounding_box = np.array(
            self.get_bounding_box_from_points()).reshape(2, -1)
        bounding_box = bounding_box[1] - bounding_box[0]

        # Define periodic mesh contraints
        for slave, master in PERIODIC_PAIRS.items():
            master_entities = self.boundary[master]
            slave_entities = self.boundary[slave]
            if len(master_entities) == 0 or len(slave_entities) == 0:
                continue

            for dim in (1, 2):
                tags_master = [j for i, j in master_entities if i == dim]
                tags_slave = [j for i, j in slave_entities if i == dim]
                gmsh.model.mesh.set_periodic(
                    dim, tags_slave, tags_master,
                    translation(SLAVE_GROUPS_OFFSET_AXES[slave], bounding_box))

        self._periodic_mesh = True

    def build_node_sets(self):
        """
        Build node sets based on the physical groups from Gmsh.

        The node sets are needed for the correct application of boundary
        conditions.

        Returns
        -------
        None.

        """

        def get_periodic_nodes(slave):
            """Get periodically connected nodes based on a slave group."""

            # Get the corresponding master group
            master = PERIODIC_PAIRS.get(slave)
            if master is None:
                slave = master

            # Compute the sorting order for finding matching nodes based on
            # the coordinates
            group_to_coords = self.boundary.get_group_coordinate_identifiers()
            order = np.where(np.isnan(group_to_coords[master]))[0][::-1]
            order = order if len(order) > 0 else [0]

            # Find all coordinates on the respective groups
            node_tags = {}
            for group in [master, slave]:
                nodes = [gmsh.model.mesh.get_nodes(dim, tag, True, False)
                         for (dim, tag) in self.boundary[group]]
                group_node_tags = [tag for n in nodes for tag in n[0]]
                group_node_coords = [coord for n in nodes for coord in n[1]]
                group_node_coords = np.reshape(group_node_coords, (-1, 3))

                # Sort the nodes according to their spatial positions
                sorter = np.lexsort([group_node_coords[:, i] for i in order])
                group_node_tags = np.array(group_node_tags, dtype=int)[sorter]

                # Make sure nodes are unique without breaking the order
                _, indices = np.unique(group_node_tags, return_index=True)
                node_tags[group] = group_node_tags[np.sort(indices)]

            master_tags = node_tags[master]
            slave_tags = node_tags[slave]
            return master_tags, slave_tags

        def remove_child_nodes(parent, children):
            """Remove nodes on the boundary of an entity from its node set."""

            for child in children:
                self.node_sets[parent] = np.setdiff1d(self.node_sets[parent],
                                                      self.node_sets[child],
                                                      assume_unique=True)
                remove_child_nodes(parent, BOUNDARY_CONNECTIVITY[child])

        # Activate the model in Gmsh
        self.activate()

        # Return if there are no nodes in the mesh
        if self.mesh_dimension < 0:
            raise RuntimeError("There are no nodes in the model.")

        # Fill up the node set for each boundary group
        for slave, master in PERIODIC_PAIRS.items():
            master_nodes, slave_nodes = get_periodic_nodes(slave)
            self.node_sets[master] = master_nodes
            self.node_sets[slave] = slave_nodes

        # Remove nodes being part of the boundary of a group
        for group in self.boundary:
            remove_child_nodes(group, BOUNDARY_CONNECTIVITY[group])

        # Add sets collecting nodes for points, edges and faces
        special_sets = [("POINT_NODES", POINT_GROUPS),
                        ("EDGE_NODES", EDGE_GROUPS),
                        ("FACE_NODES", FACE_GROUPS),
                        ("BOUNDARY_NODES", (POINT_GROUPS + EDGE_GROUPS +
                                            FACE_GROUPS))]
        for set_name, set_items in special_sets:
            self._add_set_from_sets(self.node_sets, set_name, set_items)

        # Add sets for inner nodes and all nodes nodes
        all_nodes = gmsh.model.mesh.get_nodes(gmsh.model.get_dimension(),
                                              includeBoundary=True)[0]
        self.node_sets["ALL_NODES"] = np.unique(all_nodes)
        self.node_sets["INNER_NODES"] = np.setdiff1d(
            self.node_sets["ALL_NODES"], self.node_sets["BOUNDARY_NODES"])

        # Add sets based on physical groups
        physical_groups = gmsh.model.get_physical_groups()
        for dim, tag in physical_groups:
            name = gmsh.model.get_physical_name(dim, tag)
            if len(name) == 0:
                continue
            if name not in self.boundary:
                nodes, _ = gmsh.model.mesh.get_nodes_for_physical_group(dim,
                                                                        tag)
                self.node_sets[name] = np.array(nodes, dtype=int)

    def build_element_sets(self):
        """
        Build node sets based on the physical groups from Gmsh.

        The element sets are needed to assign e.g. material to the elements.

        Returns
        -------
        None.

        """

        # Activate the model in Gmsh
        self.activate()

        # Return if there are no elements in the mesh
        if self.mesh_dimension < 0:
            raise RuntimeError("There are no elements in the model.")

        # Find all groups that need to be processed based on the mesh dimension
        relevant_groups = gmsh.model.get_physical_groups(self.mesh_dimension)

        # Initialize element sets
        element_sets = {gmsh.model.get_physical_name(*group): []
                        for group in relevant_groups}

        # Add all sets containing elements of the highest available dimension
        for group in relevant_groups:
            name = gmsh.model.get_physical_name(*group)
            if len(name) == 0:
                continue
            entities = gmsh.model.get_entities_for_physical_group(*group)
            for entity in entities:
                elements = gmsh.model.mesh.get_elements(group[0], entity)[1]
                elements = [e for el_type in elements for e in el_type]
                element_sets[name].extend(elements)

        # Make sure the elements only appear once in the sets
        for name in element_sets.keys():
            element_sets[name] = np.unique(element_sets[name])

        # Add a set for all elements in the model
        all_elements = gmsh.model.mesh.get_elements(self.mesh_dimension)[1]
        all_elements = np.array([e for el_type in all_elements
                                 for e in el_type])
        element_sets["ALL_ELEMENTS"] = all_elements

        # Assign sets to the object property and make sure elements are
        # integers
        for name in element_sets.keys():
            self.element_sets[name] = np.array(element_sets[name], dtype=int)

    def _add_set_from_sets(self, set_dictionary, set_name, set_items):
        """Add node or element sets combining multiple existing sets."""

        new_set = []
        for i in set_items:
            if set_dictionary.get(i) is not None:
                new_set.extend(set_dictionary.get(i))
        set_dictionary[set_name] = np.array(new_set, dtype=int)

    def write_abaqus_input(self, el_type_map=None):
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

        # Activate the model in Gmsh
        self.activate()

        # Load default element type map
        if el_type_map is None:
            el_type_map = ABAQUS_ELEMENT_TYPES

        # Get node information
        node_tags, node_coords, _ = gmsh.model.mesh.get_nodes(
            gmsh.model.get_dimension(), includeBoundary=True)
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

        return path

    def write_matlab_input(self):
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

        # Activate the model in Gmsh
        self.activate()

        if not self.model_dimension == 3 and self.mesh_dimension == 1:
            raise RuntimeError("Matlab support is only available for 3d "
                               "lattice unit cells.")

        # Compute RVE size
        bounding_box = np.array(
            self.get_bounding_box_from_points).reshape(2, -1)
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

        return path


class Boundary(UserDict):
    """
    The boundary information of a rectangular or box-shaped Gmsh RVE.

    Boundary inherits from collections.UserDict such that its instances offer
    the same functionalities as regular dictionaries.

    Parameters
    ----------
    model_name : str
        The name of the underlying Gmsh model.
    *args
        Arguments that will be passed to the UserDict constructor.
    sampling_density : int, optional
        The number of sampling points per spatial direction that will be
        generated to check for periodicity and validity of the model. The
        default is 10.
    periodic_default : bool
        If True, assume a periodic boundary formulation in case the automatic
        periodicity check fails. Only set this to True if you really know what
        you are doing. The default is False.
    **kwargs
        Keyword arguments that will be passed to the UserDict constructor.

    Attributes
    ----------
    periodic : bool
        A flag to mark whether the boundary is periodic.
    model_name : str
        The name of the underlying Gmsh model.
    sampling_density : bool
        The number of sampling points per spatial direction that will be
        generated to check for periodicity and validity of the model.
    data : dict[str, list[tuple[int, int]]]
        The underlying dictionary mapping the group identifier key to the
        Gmsh dimtags being part of that group. In case of periodic boundaries
        the connected dimtags have the same index in the lists.

    Example
    -------
    >>> import gmsh
    >>> from abqhom.gmshmodel import Boundary
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
    >>> # Create the boundary information
    >>> model_boundary = Boundary("SQUARE_LATTICE")
    >>> bounding_box = model_boundary.get_bounding_box_from_points()

    """

    def __init__(self, model_name, *args, sampling_density=10,
                 periodic_default=False, **kwargs):
        if not gmsh.is_initialized():
            gmsh.initialize()
        if model_name not in gmsh.model.list():
            raise ValueError("The Gmsh model {} does not exist.".format(
                model_name))

        super().__init__(*args, **kwargs)
        self._model_name = model_name
        self._sampling_density = sampling_density

        (x_min, y_min, z_min,
         x_max, y_max, z_max) = self.get_bounding_box_from_points()
        tol = gmsh.option.get_number("Geometry.Tolerance")
        model_dim = 3 - np.isclose([x_max - x_min,
                                    y_max - y_min,
                                    z_max - z_min], 0.0, atol=tol).sum()

        # Build the boundary automatically
        if len(args) == 0 and len(kwargs) == 0:
            try:
                self._build_boundary(model_dim)
            except Exception:
                raise RuntimeError(
                    "Automatic boundary generation not possible. Please "
                    "provide the boundary manually.")

        # Rearrange boundary map in case of a periodic geometry
        try:
            self._check_for_periodicity(model_dim)
        except Exception:
            self._periodic = periodic_default

        # Add missing groups as empty groups
        for key in POINT_GROUPS + EDGE_GROUPS + FACE_GROUPS:
            if key not in self:
                self[key] = []

    @property
    def periodic(self):
        """A flag to mark whether the boundary is periodic."""
        return self._periodic

    @property
    def model_name(self):
        """The name of the model the boundary is referring to."""
        return self._model_name

    @property
    def sampling_density(self):
        """The number of sample points that will be generated on entities."""
        return self._sampling_density

    @sampling_density.setter
    def sampling_density(self, n):
        """The number of sample points that will be generated on entities."""
        if n < 1:
            raise ValueError("The sampling density must be greater than zero.")
        self._sampling_density = int(n)

    def __setitem__(self, key, value):
        if (not isinstance(key, str) or
                key.upper() not in set(POINT_GROUPS +
                                       EDGE_GROUPS +
                                       FACE_GROUPS)):
            raise KeyError("Key must be a group identifier string.")
        super().__setitem__(key.upper(), value)

    def __getitem__(self, key):
        if (not isinstance(key, str) or
                key.upper() not in set(POINT_GROUPS +
                                       EDGE_GROUPS +
                                       FACE_GROUPS)):
            raise KeyError("Key must be a group identifier string.")
        return super().__getitem__(key.upper())

    def get_bounding_box_from_points(self):
        """
        Compute the bounding box of the geometric CAD model from points.

        Returns
        -------
        x_min : float
            The minimum x-value of the bounding box.
        y_min : float
            The minimum y-value of the bounding box.
        z_min : float
            The minimum z-value of the bounding box.
        x_max : float
            The maximum x-value of the bounding box.
        y_max : float
            The maximum y-value of the bounding box.
        z_max : float
            The maximum z-value of the bounding box.

        """

        # This method is mainly needed since the bounding box Gmsh computes
        # by using gmsh.model.get_bounding_box() delivers highly inaccurate
        # results, especially when using the OCC kernel for geometry design

        gmsh.model.set_current(self.model_name)

        # Get all point entities
        points = gmsh.model.get_entities(0)
        if len(points) == 0:
            return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

        # Compare coordinates of points to find minimum and maximum
        min_coords = gmsh.model.get_value(points[0][0], points[0][1], [])
        max_coords = gmsh.model.get_value(points[0][0], points[0][1], [])
        min_coords, max_coords = np.array(min_coords), np.array(max_coords)
        for p in points:
            coords = np.array(gmsh.model.get_value(p[0], p[1], []))
            min_coords[min_coords > coords] = coords[min_coords > coords]
            max_coords[max_coords < coords] = coords[max_coords < coords]

        # Extract minimum and maximum from arrays
        x_min, y_min, z_min = min_coords
        x_max, y_max, z_max = max_coords

        return x_min, y_min, z_min, x_max, y_max, z_max

    def _build_boundary(self, model_dim):
        """Build the boundary information for the requested Gmsh model."""

        # Activate the model
        gmsh.model.set_current(self.model_name)

        bounding_box = np.array(self.get_bounding_box_from_points())
        bounding_box = bounding_box.reshape(2, -1)
        delta_x, delta_y, delta_z = bounding_box[1] - bounding_box[0]
        dim = 3 - np.isclose([delta_x, delta_y, delta_z], 0.0, atol=1e-6).sum()

        # Get the model bounding box
        tol = gmsh.option.get_number("Geometry.Tolerance")
        group_to_coords = self.get_group_coordinate_identifiers()
        (x_min, y_min, z_min,
         x_max, y_max, z_max) = self.get_bounding_box_from_points()

        # Initialize boundary map
        for group_name in group_to_coords.keys():
            self[group_name] = []

        # Start looping over all model entities
        entities = gmsh.model.get_entities()
        for dim, tag in entities:
            if dim == 3:
                continue

            # Generate sample points
            coords = self._generate_sample_points_on_entity(
                dim, tag, self.sampling_density)

            # Check if the entity is entirely inside the bounding box
            x, y, z = coords[:, 0], coords[:, 1], coords[:, 2]
            for k, k_min, k_max in zip([x, y, z],
                                       [x_min, y_min, z_min],
                                       [x_max, y_max, z_max]):
                min_check = all((k > k_min) | np.isclose(k, k_min, atol=tol))
                max_check = all((k < k_max) | np.isclose(k, k_max, atol=tol))
                if (not min_check) or (not max_check):
                    raise RuntimeError(
                        "Entity ({}, {}) is outside of the computed bounding "
                        "box".format(dim, tag))

            # Check if entities are part of a boundary group
            for group_name, identifiers in group_to_coords.items():
                group_type = group_name.split("_")[0]
                if group_type in ["FACE", "EDGE", "POINT"][:(3 - model_dim)]:
                    continue

                mask = np.invert(np.isnan(identifiers))
                if np.allclose(coords[:, mask], identifiers[mask], atol=tol):
                    self[group_name].append((dim, tag))
                    break

        # Cleanup of boundary groups
        for group_name, group_dimtags in self.items():
            to_remove = set()
            for dimtag in group_dimtags:
                upward, downward = gmsh.model.get_adjacencies(*dimtag)

                # Points such as spline control points will be removed
                if len(upward) == len(downward) == 0:
                    to_remove.add(dimtag)

                # Entities that are on the boundary of higher order entities
                # in the same group will be removed as well
                downward = {(dimtag[0] - 1, tag) for tag in downward}
                to_remove.update(downward & set(group_dimtags))

            self[group_name] = list(set(group_dimtags) - to_remove)

    def get_group_coordinate_identifiers(self):
        """Get all spatial coordinates being part of a point, edge or face."""
        # Create coordinate identifiers for boundary groups
        (x_min, y_min, z_min,
         x_max, y_max, z_max) = self.get_bounding_box_from_points()
        group_identifiers = np.array([[x_min, y_min, z_min],     # POINT_1
                                      [x_max, y_min, z_min],     # POINT_2
                                      [x_max, y_max, z_min],     # POINT_3
                                      [x_min, y_max, z_min],     # POINT_4
                                      [x_min, y_min, z_max],     # POINT_5
                                      [x_max, y_min, z_max],     # POINT_6
                                      [x_max, y_max, z_max],     # POINT_7
                                      [x_min, y_max, z_max],     # POINT_8
                                      [np.nan, y_min, z_min],    # EDGE_1
                                      [x_max, np.nan, z_min],    # EDGE_2
                                      [np.nan, y_max, z_min],    # EDGE_3
                                      [x_min, np.nan, z_min],    # EDGE_4
                                      [np.nan, y_min, z_max],    # EDGE_5
                                      [x_max, np.nan, z_max],    # EDGE_6
                                      [np.nan, y_max, z_max],    # EDGE_7
                                      [x_min, np.nan, z_max],    # EDGE_8
                                      [x_min, y_min, np.nan],    # EDGE_9
                                      [x_max, y_min, np.nan],    # EDGE_10
                                      [x_max, y_max, np.nan],    # EDGE_11
                                      [x_min, y_max, np.nan],    # EDGE_12
                                      [np.nan, np.nan, z_min],   # FACE_1
                                      [np.nan, np.nan, z_max],   # FACE_2
                                      [np.nan, y_min, np.nan],   # FACE_3
                                      [np.nan, y_max, np.nan],   # FACE_4
                                      [x_min, np.nan, np.nan],   # FACE_5
                                      [x_max, np.nan, np.nan]])  # FACE_6

        # Assign the identifiers to the groups names with ascending order
        group_to_coords = {}
        for i, j in enumerate(POINT_GROUPS + EDGE_GROUPS + FACE_GROUPS):
            group_to_coords[j] = group_identifiers[i]

        return group_to_coords

    def _generate_sample_points_on_entity(self, dim, tag, n):
        """Create n sample points on a Gmsh entitiy."""
        (bounds_min,
         bounds_max) = gmsh.model.get_parametrization_bounds(dim, tag)
        sample_points = [np.linspace(i, j, n)
                         for i, j in zip(bounds_min, bounds_max)]
        sample_points = np.meshgrid(*sample_points)

        # Get the physical coordinates from the parametric space
        assert len(sample_points) in {0, 1, 2}
        if len(sample_points) == 0:  # Point entities
            sample_points = []
        elif len(sample_points) == 1:  # Line entities
            sample_points = sample_points[0]
        else:  # Face entities
            u, v = sample_points
            sample_points = np.vstack((u.flatten(), v.flatten())).T.flatten()
        coords = gmsh.model.get_value(dim, tag, sample_points)
        coords = coords.reshape(-1, 3)

        # Filter sample points of trimmed entities to make sure to only
        # obtain points inside the trimmed domain
        if dim > 0:
            inside = np.array([gmsh.model.is_inside(dim, tag, c)
                               for c in coords])
            coords = coords[inside.nonzero()]

        return coords

    def _check_for_periodicity(self, model_dim):
        """If possible, reorder the boundary map for periodic constraints."""

        def sort_sample_coordinates(coords, master):
            """Sort coordinate arrays according to the below dictionary."""
            group_to_coords = self.get_group_coordinate_identifiers()
            order = np.where(np.isnan(group_to_coords[master]))[0][::-1]
            order = order if len(order) > 0 else [0]
            coords = coords[np.lexsort([coords[:, i] for i in order])]
            return coords

        boundary_groups = {0: ["POINT_1"],
                           1: ["POINT_1", "POINT_2"],
                           2: ["POINT_1", "POINT_2", "POINT_3", "POINT_4",
                               "EDGE_1", "EDGE_2", "EDGE_3", "EDGE_4"],
                           3: POINT_GROUPS + EDGE_GROUPS + FACE_GROUPS}

        # Make a deepcopy of the original boundary map
        original_boundary_map = {i: j.copy() for i, j in self.data.items()}

        self._periodic = True
        tol = gmsh.option.get_number("Geometry.Tolerance")
        group_to_coords = self.get_group_coordinate_identifiers()
        for slave, master in PERIODIC_PAIRS.items():
            slave_indices = []
            master_entities = self[master]
            slave_entities = self[slave]

            # Skip irrelevant slave groups

            if slave not in boundary_groups[model_dim]:
                continue

            for master_dimtag in master_entities:
                # Generate sample points on the master entity
                coords_master = self._generate_sample_points_on_entity(
                    *master_dimtag, self.sampling_density)
                coords_master = sort_sample_coordinates(coords_master, master)

                # Check periodicity for all potential slave entities
                for slave_index, slave_dimtag in enumerate(slave_entities):
                    if master_dimtag[0] != slave_dimtag[0]:
                        continue

                    # Generate sample points on the potential slave entity
                    coords_slave = self._generate_sample_points_on_entity(
                        *slave_dimtag, self.sampling_density)
                    coords_slave = sort_sample_coordinates(coords_slave,
                                                           master)

                    # If the number of points generated is not the same the
                    # entities cannot be periodically connected
                    if len(coords_master) != len(coords_slave):
                        continue

                    # Check if the master and slave are periodically connected
                    delta_true = np.abs(coords_slave - coords_master)
                    delta_expected = np.nan_to_num(np.abs(
                        group_to_coords[slave] - group_to_coords[master]))
                    if np.allclose(delta_true, delta_expected, atol=tol):
                        slave_indices.append(slave_index)
                        break

            # Update the order of the slave entities and check periodicity
            self[slave] = [slave_entities[i] for i in slave_indices]
            if len(self[slave]) != len(master_entities):
                self._periodic = False
                self.data = original_boundary_map
                return
