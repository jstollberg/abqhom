# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-or-later

"""Gmsh RVE geometry, boundary classification, and mesh export."""

import os
from collections import OrderedDict

try:
    from collections import UserDict
except ImportError:
    from UserDict import UserDict

import numpy as np
import gmsh

from abqhom.constants import (POINT_GROUPS, EDGE_GROUPS, FACE_GROUPS,
                              PERIODIC_PAIRS, ABAQUS_ELEMENT_TYPES_2D,
                              ABAQUS_ELEMENT_TYPES_3D, BOUNDARY_CONNECTIVITY,
                              ABAQUS_CONNECTIVITY_ORDER,
                              ELEMENTS_NOT_SUPPORTED,
                              ELEMENTS_WARPING_DOFS)


class GmshRVE(object):
    """
    Gmsh representation of a structural or continuum RVE.

    The reference cell must be rectangular in the global X1-X2 plane with
    X3 = 0 in 2D or cuboidal and aligned with the global coordinate axes in
    3D. Boundary entities and their periodic primary--secondary relations are
    detected automatically unless a :class:`Boundary` is supplied. Physical
    groups define the continuum, truss, beam, and shell constituents exported
    to Abaqus; those constituents may have arbitrary orientations within or
    on the cell boundary.

    Parameters
    ----------
    model_name : str
        The name of the underlying Gmsh model.
    boundary : Boundary, optional
        Boundary-entity map. If None, it is generated automatically.

    Attributes
    ----------
    model_name : str
        The name of the underlying Gmsh model.
    boundary : Boundary
        The boundary definition of the RVE.
    periodic_boundary : bool
        A flag to mark whether the model is periodic.
    periodic_mesh : bool
        A flag indicating whether periodic mesh constraints were prescribed
        or opposite boundary node sets were detected as matching.
    node_sets : dict[str, ndarray]
        A map from set names to arrays of Gmsh node tags for user-defined
        physical groups and further default node sets.
    element_sets : dict[str, ndarray]
        A map from set names to arrays of Gmsh element tags for user-defined
        physical groups and further default element sets.
    model_dimension : int
        The geometric dimensionality of the model.
    mesh_dimension : int
        The dimensionality of the mesh of the model.

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
        self._periodic_mesh = False

        # Create the boundary map
        if boundary is None:
            boundary = Boundary(model_name)
        elif boundary.model_name != model_name:
            raise ValueError(
                "The supplied boundary belongs to Gmsh model '{}', not "
                "'{}'.".format(boundary.model_name, model_name))
        self._boundary = boundary

        # Check if the model has entities of at least dimension 1
        self.activate()
        if len(gmsh.model.get_entities(1)) == 0:
            raise RuntimeError(
                "The model must contain entities of at least dimension 1.")

        self._validate_model_geometry()

        # Make sure the RVE actually has a detectable boundary
        all_boundary_groups = POINT_GROUPS + EDGE_GROUPS + FACE_GROUPS
        if all(len(self.boundary[g]) == 0 for g in all_boundary_groups):
            raise RuntimeError(
                "No boundary entities were detected on the RVE. If the "
                "automatic detection fails, provide the boundary manually "
                "through the optional 'boundary' argument.")

        # Add physical boundary groups to the model
        self._add_boundary_groups()

    @property
    def model_name(self):
        """
        The name of the model.

        Returns
        -------
        str
            The model name.

        """

        return self.boundary.model_name

    @property
    def boundary(self):
        """
        The boundary definition of the RVE.

        Returns
        -------
        Boundary
            The boundary of the RVE.

        """

        return self._boundary

    @property
    def periodic_boundary(self):
        """
        A flag to mark whether the model is periodic.

        Returns
        -------
        bool
            True if the model is periodic, False otherwise.

        """

        return self._boundary.periodic

    @property
    def periodic_mesh(self):
        """
        Whether the mesh is configured or detected as periodically matching.

        Returns
        -------
        bool
            True if periodic mesh constraints were prescribed or opposite
            boundary nodes were paired one to one, False otherwise.

        """

        return self._periodic_mesh

    @property
    def node_sets(self):
        """
        Node-tag sets, e.g. for the nodes of a boundary group.

        Returns
        -------
        dict[str, ndarray]
            Keys are set names and values are arrays of Gmsh node tags.

        """

        return self._node_sets

    @property
    def element_sets(self):
        """
        Element-tag sets, e.g. for the elements of a physical group.

        Returns
        -------
        dict[str, ndarray]
            Keys are set names and values are arrays of Gmsh element tags.

        """

        return self._element_sets

    @property
    def model_dimension(self):
        """
        The geometric dimensionality of the RVE model.

        Returns
        -------
        dim : int
            The model dimension. If no entities were created -1 is
            returned.

        """

        self.activate()

        if len(gmsh.model.get_entities()) == 0:
            return -1

        return self.boundary.model_dimension

    @property
    def mesh_dimension(self):
        """
        The dimensionality of the mesh of the RVE model.

        Returns
        -------
        dim : int
            The mesh dimension. If no mesh exists, -1 is returned.

        """

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
        """Set the RVE model as the current model in Gmsh."""

        gmsh.model.set_current(self.model_name)

    def show(self):
        """Open the RVE model in the Gmsh graphical user interface."""

        self.activate()
        gmsh.fltk.run()

    def get_bounding_box_from_points(self):
        """
        Return the RVE bounding box from CAD points or mesh nodes.

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

    def _validate_model_geometry(self):
        """Validate the supported RVE dimension and coordinate plane."""

        model_dimension = self.model_dimension
        if model_dimension not in {2, 3}:
            raise RuntimeError(
                "The RVE must be two- or three-dimensional.")

        # Abaqus two-dimensional parts use the global x-y plane. Requiring
        # the Gmsh model to use the same plane avoids an implicit coordinate
        # transformation during export.
        if model_dimension == 2:
            bounding_box = np.asarray(
                self.get_bounding_box_from_points(),
                dtype=np.float64).reshape(2, 3)
            tolerance = gmsh.option.get_number("Geometry.Tolerance")
            if not np.allclose(
                    bounding_box[:, 2], 0.0, atol=tolerance, rtol=0.0):
                raise RuntimeError(
                    "A two-dimensional RVE must lie in the global x-y "
                    "plane (z = 0).")

    def get_nodes(self):
        """
        Return constituent mesh nodes in the original Gmsh coordinates.

        Automatically generated boundary physical groups are excluded.

        Returns
        -------
        node_tags : ndarray
            Node tags.
        node_coords : ndarray
            Node coordinates.

        """

        self.activate()

        # Get nodes from physical groups
        node_tags = []
        node_coords = []
        for group in gmsh.model.get_physical_groups():
            if gmsh.model.get_physical_name(*group) in (POINT_GROUPS
                                                        + EDGE_GROUPS
                                                        + FACE_GROUPS):
                continue

            entities = gmsh.model.get_entities_for_physical_group(*group)
            for tag in entities:
                (entity_node_tags,
                 entity_node_coords,
                 _) = gmsh.model.mesh.get_nodes(group[0], tag,
                                                includeBoundary=True)
                node_tags.append(entity_node_tags)
                node_coords.append(entity_node_coords)

        # Combine nodes and coordinates from different sources
        if len(node_tags) == 0:
            raise RuntimeError(
                "No constituent physical groups containing mesh nodes were "
                "found.")
        node_tags = np.concatenate(node_tags)
        node_coords = np.concatenate(node_coords)
        node_coords = np.array(node_coords).reshape(-1, 3)

        # Make sure to not have duplicate nodes from physical groups
        _, unique_idx = np.unique(node_tags, return_index=True)
        node_tags = node_tags[unique_idx]
        node_coords = node_coords[unique_idx]

        return node_tags, node_coords

    def get_elements(self):
        """
        Return all mesh elements belonging to constituent physical groups.

        Automatically generated boundary physical groups are excluded.

        Returns
        -------
        element_tags : dict[int, ndarray]
            The relevant element tags. Keys are the Gmsh element types and
            values are the element tags of the respective element type.
        element_connectivity : dict[int, ndarray]
            The relevant element connectivity lists. Keys are the Gmsh element
            types and values are the connectivity lists.

        """

        self.activate()

        element_tags = {}
        element_connectivity = {}

        # Get element information from physical groups
        for group in gmsh.model.get_physical_groups():
            if gmsh.model.get_physical_name(*group) in (POINT_GROUPS
                                                        + EDGE_GROUPS
                                                        + FACE_GROUPS):
                continue

            entities = gmsh.model.get_entities_for_physical_group(*group)
            for tag in entities:
                (entity_el_types,
                 entity_el_tags,
                 entity_conn) = gmsh.model.mesh.get_elements(group[0], tag)

                for i, el_type in enumerate(entity_el_types):
                    element_tags[el_type] = element_tags.get(el_type, [])
                    element_connectivity[el_type] = element_connectivity.get(
                        el_type, [])
                    element_tags[el_type].append(entity_el_tags[i])
                    element_connectivity[el_type].append(entity_conn[i])

        if len(element_tags) == 0:
            raise RuntimeError(
                "No constituent physical groups containing mesh elements "
                "were found.")

        # Rearrange connectivity list
        for el_type in element_tags:
            element_tags[el_type] = np.concatenate(element_tags[el_type])
            conn = np.concatenate(element_connectivity[el_type])

            nen = gmsh.model.mesh.getElementProperties(el_type)[3]
            conn = conn.reshape(-1, nen)
            if el_type not in ABAQUS_CONNECTIVITY_ORDER:
                raise NotImplementedError(
                    "Gmsh element type {} has no supported Abaqus "
                    "connectivity mapping.".format(el_type))
            conn = conn[:, ABAQUS_CONNECTIVITY_ORDER[el_type]]
            element_connectivity[el_type] = conn

        # Remove duplicate elements
        for el_type in element_tags:
            tags = element_tags[el_type]
            conn = element_connectivity[el_type]

            _, unique_idx = np.unique(element_tags[el_type], return_index=True)
            element_tags[el_type] = tags[unique_idx]
            element_connectivity[el_type] = conn[unique_idx]

        return element_tags, element_connectivity

    def _add_boundary_groups(self):
        """Add physical groups for boundary entities to the Gmsh model."""

        self.activate()

        # Check for name conflicts
        boundary_groups = POINT_GROUPS + EDGE_GROUPS + FACE_GROUPS
        for dimtag in gmsh.model.get_physical_groups():
            name = gmsh.model.get_physical_name(*dimtag)
            if name.upper() in boundary_groups:
                raise RuntimeError(
                    "Physical group {} conflicts with the automatically "
                    "generated physical boundary groups. Do not use any of "
                    "the following names: {}".format(name, boundary_groups))

        # Add physical groups to model
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
        element_order : int, optional
            The order of the generated elements. Possible values are 1
            and 2. The default is 1.
        algorithm_2d : int, optional
            The 2D meshing algorithm to use. The available algorithms
            are listed in the Gmsh documentation. The default is 6.
        algorithm_3d : int, optional
            The 3D meshing algorithm to use. The available algorithms
            are listed in the Gmsh documentation. The default is 1.
        element_size : float, optional
            The characteristic element length. If this is None Gmsh will
            set the element size automatically. The default is None.

        """

        self.activate()

        # Set element order
        if element_order not in {1, 2}:
            raise ValueError("Only first-order and second-order elements "
                             "are supported.")
        gmsh.option.set_number("Mesh.ElementOrder", element_order)
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
        Impose matching meshes on corresponding periodic boundary entities.

        Secondary entities are mapped to their primary partners by the
        appropriate primitive-cell translation.

        """

        def translation(offset):
            """Return the affine matrix for a translation vector."""

            transformation = np.array([1.0, 0.0, 0.0, 0.0,
                                       0.0, 1.0, 0.0, 0.0,
                                       0.0, 0.0, 1.0, 0.0,
                                       0.0, 0.0, 0.0, 1.0])
            for axis, value in enumerate(offset):
                transformation[axis*4 + 3] = value
            return transformation

        if not self.periodic_boundary:
            raise RuntimeError(
                "The model does not have a periodic boundary formulation.")

        group_to_coords = self.boundary.get_group_coordinate_identifiers()

        # Define periodic mesh constraints
        for secondary, primary in PERIODIC_PAIRS.items():
            primary_entities = self.boundary[primary]
            secondary_entities = self.boundary[secondary]
            if len(primary_entities) == 0 or len(secondary_entities) == 0:
                continue

            for dim in (1, 2):
                primary_tags = [j for i, j in primary_entities if i == dim]
                secondary_tags = [j for i, j in secondary_entities
                                  if i == dim]
                offset = np.nan_to_num(
                    group_to_coords[secondary]
                    - group_to_coords[primary])
                gmsh.model.mesh.set_periodic(
                    dim, secondary_tags, primary_tags,
                    translation(offset))

        self._periodic_mesh = True

    def build_node_sets(self):
        """
        Build constituent and disjoint boundary node sets.

        Nodes on intersections of boundary entities are assigned to the
        lowest-dimensional group. Periodic primary and secondary sets are
        ordered so that nodes with the same list index form a pair.

        """

        def get_periodic_nodes(secondary):
            """Return ordered node pairs for a secondary boundary group."""

            # Get the corresponding primary group
            primary = PERIODIC_PAIRS.get(secondary)
            if primary is None:
                primary = secondary

            # Compute the sorting order for finding matching nodes based
            # on the coordinates
            group_to_coords = self.boundary.get_group_coordinate_identifiers()
            order = np.where(np.isnan(group_to_coords[primary]))[0][::-1]
            order = order if len(order) > 0 else [0]

            # Find all coordinates on the respective groups
            node_tags = {}
            for group in [primary, secondary]:
                nodes = [gmsh.model.mesh.get_nodes(dim, tag, True, False)
                         for (dim, tag) in self.boundary[group]]
                group_node_tags = [tag for n in nodes for tag in n[0]]
                group_node_coords = [coord for n in nodes for coord in n[1]]
                group_node_coords = np.reshape(group_node_coords, (-1, 3))
                group_node_coords = np.round(group_node_coords, decimals=10)

                # Sort the nodes according to their spatial positions
                sorter = np.lexsort([group_node_coords[:, i] for i in order])
                group_node_tags = np.array(group_node_tags,
                                           dtype=np.uint64)[sorter]

                # Make sure nodes are unique without breaking the order
                _, indices = np.unique(group_node_tags, return_index=True)
                node_tags[group] = group_node_tags[np.sort(indices)]

            primary_tags = node_tags[primary]
            secondary_tags = node_tags[secondary]
            return primary_tags, secondary_tags

        def remove_child_nodes(parent, children):
            """Remove nodes of all child boundary groups from a parent set."""

            for child in children:
                self.node_sets[parent] = np.setdiff1d(self.node_sets[parent],
                                                      self.node_sets[child],
                                                      assume_unique=True)
                remove_child_nodes(parent, BOUNDARY_CONNECTIVITY[child])

        # Activate the model in Gmsh
        self.activate()

        # Reject models without mesh nodes.
        if self.mesh_dimension < 0:
            raise RuntimeError("There are no nodes in the model.")

        # Fill up the node set for each boundary group
        for secondary, primary in PERIODIC_PAIRS.items():
            paired_nodes = get_periodic_nodes(secondary)
            self.node_sets[primary] = paired_nodes[0]
            self.node_sets[secondary] = paired_nodes[1]

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
            self._new_set_from_union(self.node_sets, set_name, *set_items)

        # Add sets for inner nodes and all nodes
        all_nodes, _ = self.get_nodes()
        self.node_sets["ALL_NODES"] = np.unique(all_nodes).astype(np.uint64)
        self.node_sets["INNER_NODES"] = np.setdiff1d(
            self.node_sets["ALL_NODES"], self.node_sets["BOUNDARY_NODES"])

        # Add sets based on physical groups
        for dim, tag in gmsh.model.get_physical_groups():
            name = gmsh.model.get_physical_name(dim, tag)
            if len(name) == 0:
                continue
            if name not in self.boundary.keys():
                nodes, _ = gmsh.model.mesh.get_nodes_for_physical_group(dim,
                                                                        tag)
                self.node_sets[name] = np.array(nodes, dtype=np.uint64)

    def _get_boundary_signatures(self):
        """Return coordinate signatures of the RVE boundary groups."""

        bounding_box = np.array(
            self.get_bounding_box_from_points()).reshape(2, -1)
        bounds_min, bounds_max = bounding_box
        # Supported 2D cells occupy x-y; 3D cells use all coordinate axes.
        active_axes = np.arange(3) < self.model_dimension

        boundary_groups = POINT_GROUPS + EDGE_GROUPS + FACE_GROUPS
        if self.model_dimension == 2:
            boundary_groups = POINT_GROUPS[:4] + EDGE_GROUPS[:4]

        group_to_coords = self.boundary.get_group_coordinate_identifiers()
        group_signature, signature_to_group = {}, {}
        for group in boundary_groups:
            signature = []
            for axis, value in enumerate(group_to_coords[group]):
                if (not active_axes[axis]) or np.isnan(value):
                    signature.append(0)
                elif (abs(value - bounds_min[axis])
                      < abs(value - bounds_max[axis])):
                    signature.append(-1)
                else:
                    signature.append(1)
            group_signature[group] = tuple(signature)
            signature_to_group[tuple(signature)] = group

        return group_signature, signature_to_group, active_axes

    def _get_boundary_element_sets(self, element_tags,
                                   element_connectivity, tol=None):
        """Classify elements whose complete nodal support is on a boundary."""

        bounding_box = np.array(
            self.get_bounding_box_from_points()).reshape(2, -1)
        bounds_min, bounds_max = bounding_box
        if tol is None:
            tol = 1e-6*float(np.max(bounds_max - bounds_min)) + 1e-12

        _, signature_to_group, active_axes = self._get_boundary_signatures()

        node_tags, node_coords = self.get_nodes()
        node_indices = {int(tag): i for i, tag in enumerate(node_tags)}
        boundary_groups = POINT_GROUPS + EDGE_GROUPS + FACE_GROUPS
        element_sets = {group: [] for group in boundary_groups}

        for el_type, tags in element_tags.items():
            connectivity = element_connectivity[el_type]
            for tag, nodes in zip(tags, connectivity):
                indices = [node_indices[int(node)] for node in nodes]
                coords = node_coords[indices]

                signature = []
                for axis in range(3):
                    if not active_axes[axis]:
                        signature.append(0)
                    elif np.all(np.abs(coords[:, axis]
                                       - bounds_min[axis]) < tol):
                        signature.append(-1)
                    elif np.all(np.abs(coords[:, axis]
                                       - bounds_max[axis]) < tol):
                        signature.append(1)
                    else:
                        signature.append(0)

                group = signature_to_group.get(tuple(signature))
                if group is not None:
                    element_sets[group].append(tag)

        for group, tags in element_sets.items():
            element_sets[group] = np.array(tags, dtype=np.uint64)

        return element_sets

    def build_element_sets(self):
        """
        Build constituent, boundary, and reduction element sets.

        In addition to the physical groups, the method creates
        ``BOUNDARY_ELEMENTS`` and ``INNER_ELEMENTS``. Elements whose complete
        nodal support belongs to a secondary periodic entity are collected in
        ``SECONDARY_ELEMENTS`` for the ownership-consistent reduction.

        """

        # Activate the model in Gmsh
        self.activate()

        # Reject models without mesh elements.
        if self.mesh_dimension < 0:
            raise RuntimeError("There are no elements in the model.")

        # Add element sets based on the physical groups
        for group in gmsh.model.get_physical_groups():
            name = gmsh.model.get_physical_name(*group)
            if (len(name) == 0) or (name in (POINT_GROUPS + EDGE_GROUPS
                                             + FACE_GROUPS)):
                continue

            group_elements = []
            entities = gmsh.model.get_entities_for_physical_group(*group)
            for entity in entities:
                elements = gmsh.model.mesh.get_elements(group[0], entity)[1]
                elements = [e for el_type in elements for e in el_type]
                group_elements.extend(elements)
            self.element_sets[name] = (
                np.unique(group_elements).astype(np.uint64))

        # Add a set for all elements in the model
        element_tags, element_connectivity = self.get_elements()
        all_elements = [e for t in element_tags.values() for e in t]
        self.element_sets["ALL_ELEMENTS"] = np.array(all_elements,
                                                     dtype=np.uint64)

        # Add element sets for the RVE boundary groups
        boundary_sets = self._get_boundary_element_sets(
            element_tags, element_connectivity)
        self.element_sets.update(boundary_sets)

        # Add sets collecting elements for points, edges and faces
        special_sets = [("POINT_ELEMENTS", POINT_GROUPS),
                        ("EDGE_ELEMENTS", EDGE_GROUPS),
                        ("FACE_ELEMENTS", FACE_GROUPS),
                        ("BOUNDARY_ELEMENTS", (POINT_GROUPS + EDGE_GROUPS
                                               + FACE_GROUPS)),
                        ("SECONDARY_ELEMENTS", PERIODIC_PAIRS.keys())]
        for set_name, set_items in special_sets:
            self._new_set_from_union(self.element_sets, set_name, *set_items)

        # Add a set for all elements not located entirely on the boundary
        self.element_sets["INNER_ELEMENTS"] = np.setdiff1d(
            self.element_sets["ALL_ELEMENTS"],
            self.element_sets["BOUNDARY_ELEMENTS"])

    @classmethod
    def read_from_mesh(cls, model_name):
        """
        Build an RVE from an existing two- or three-dimensional mesh.

        Boundary nodes are classified from an axis-aligned nodal bounding
        box. A 2D mesh must lie in the global X1-X2 plane with X3 = 0.
        Opposite groups are periodic when their free coordinates and
        multiplicities match one to one. Named physical groups must identify
        the constituent elements.

        Parameters
        ----------
        model_name : str
            The name of the Gmsh model holding the mesh.

        Returns
        -------
        rve : GmshRVE
            The RVE built from the existing mesh.

        """

        cell = cls.__new__(cls)
        cell._node_sets = {}
        cell._element_sets = {}
        cell._periodic_mesh = False
        empty = {g: [] for g in (POINT_GROUPS + EDGE_GROUPS + FACE_GROUPS)}
        cell._boundary = Boundary(model_name, initial_data=empty,
                                  check_periodic=False)
        cell._validate_model_geometry()
        cell._build_node_sets_from_mesh()
        cell.build_element_sets()
        return cell

    @staticmethod
    def _boundary_node_members(coords, bounds_min, bounds_max, tol,
                               signature_to_group, active_axes):
        """Classify nodal indices by their bounding-box signature."""

        flags = np.zeros((len(coords), 3), dtype=int)
        for axis in np.where(active_axes)[0]:
            flags[np.abs(coords[:, axis] - bounds_min[axis]) < tol,
                  axis] = -1
            flags[np.abs(coords[:, axis] - bounds_max[axis]) < tol,
                  axis] = 1

        members = {}
        for index, signature in enumerate(map(tuple, flags)):
            if signature == (0, 0, 0):
                continue
            group = signature_to_group.get(signature)
            if group is not None:
                members.setdefault(group, []).append(index)
        return members

    def _order_boundary_node_sets(self, tags, coords, members,
                                  group_signature, active_axes):
        """Order paired boundary nodes and return their free coordinates."""

        group_coordinates = {}
        for group, indices in members.items():
            indices = np.array(indices, dtype=int)
            free_axes = [axis for axis in range(3)
                         if (active_axes[axis]
                             and group_signature[group][axis] == 0)]
            if len(free_axes) > 0:
                keys = [np.round(coords[indices, axis], 8)
                        for axis in free_axes]
                indices = indices[np.lexsort(keys)]

            self.node_sets[group] = tags[indices]
            group_coordinates[group] = tuple(
                tuple(round(float(coords[index, axis]), 8)
                      for axis in free_axes)
                for index in indices)
        return group_coordinates

    @staticmethod
    def _periodic_node_sets_match(group_coordinates, group_signature):
        """Return whether all relevant periodic node sets match one to one."""

        for secondary, primary in PERIODIC_PAIRS.items():
            if secondary not in group_signature:
                continue
            if (group_coordinates.get(primary, ())
                    != group_coordinates.get(secondary, ())):
                return False
        return True

    def _add_prebuilt_mesh_node_sets(self, tags):
        """Add aggregate and constituent node sets for a pre-built mesh."""

        special = [("POINT_NODES", POINT_GROUPS),
                   ("EDGE_NODES", EDGE_GROUPS),
                   ("FACE_NODES", FACE_GROUPS),
                   ("BOUNDARY_NODES",
                    POINT_GROUPS + EDGE_GROUPS + FACE_GROUPS)]
        for set_name, groups in special:
            self._new_set_from_union(self.node_sets, set_name, *groups)
        self.node_sets["ALL_NODES"] = np.unique(tags).astype(np.uint64)
        self.node_sets["INNER_NODES"] = np.setdiff1d(
            self.node_sets["ALL_NODES"], self.node_sets["BOUNDARY_NODES"])

        boundary_groups = POINT_GROUPS + EDGE_GROUPS + FACE_GROUPS
        for dim, tag in gmsh.model.get_physical_groups():
            name = gmsh.model.get_physical_name(dim, tag)
            if len(name) == 0 or name in boundary_groups:
                continue
            nodes, _ = gmsh.model.mesh.get_nodes_for_physical_group(dim, tag)
            self.node_sets[name] = np.array(nodes, dtype=np.uint64)

    def _build_node_sets_from_mesh(self, tol=None):
        """Classify and pair boundary nodes of a pre-built mesh."""

        self.activate()
        tags, coords, _ = gmsh.model.mesh.get_nodes()
        tags = np.array(tags, dtype=np.uint64)
        coords = np.array(coords, dtype=float).reshape(-1, 3)

        bbox = np.array(self.get_bounding_box_from_points()).reshape(2, 3)
        lo, hi = bbox[0], bbox[1]
        if tol is None:
            tol = 1e-6*float(np.max(hi - lo))

        (group_signature,
         signature_to_group,
         active_axes) = self._get_boundary_signatures()
        members = self._boundary_node_members(
            coords, lo, hi, tol, signature_to_group, active_axes)
        group_coordinates = self._order_boundary_node_sets(
            tags, coords, members, group_signature, active_axes)
        self._periodic_mesh = self._periodic_node_sets_match(
            group_coordinates, group_signature)
        self._add_prebuilt_mesh_node_sets(tags)

    @staticmethod
    def _new_set_from_union(set_dictionary, set_name, *set_items):
        """Add the unordered union of existing node or element sets."""

        union = set()
        for i in set_items:
            union.update(set_dictionary.get(i, []))

        set_dictionary[set_name] = np.array(list(union), dtype=np.uint64)

    def write_abaqus_input(self, custom_types=None):
        """
        Export the Gmsh mesh to Abaqus.

        Node and element sets will not be written to the input file and
        need to be added to the model subsequently directly within
        Abaqus.

        Parameters
        ----------
        custom_types : mapping of int to str, optional
            Abaqus element names that override the defaults for supported
            Gmsh connectivity types. The mapping does not introduce new
            connectivity permutations. If None, the default map is used.

        Returns
        -------
        path : str
            Absolute path to the generated Abaqus input file.

        """

        # Activate the model in Gmsh
        self.activate()

        # Load default element type map
        custom_types = {} if custom_types is None else dict(custom_types)
        assert self.model_dimension in {2, 3}
        el_type_map = dict(ABAQUS_ELEMENT_TYPES_3D)
        if self.model_dimension == 2:
            el_type_map = dict(ABAQUS_ELEMENT_TYPES_2D)
        el_type_map.update(custom_types)

        # Get node information
        node_tags, node_coords = self.get_nodes()
        element_tags, conn = self.get_elements()

        # Reject element types the homogenization does not support
        used_types = set()
        for el_type in element_tags:
            if el_type not in el_type_map:
                raise NotImplementedError(
                    "No Abaqus element name is configured for Gmsh element "
                    "type {}. Supply one through 'custom_types'."
                    .format(el_type))
            used_types.add(el_type_map[el_type])
        warping = used_types & ELEMENTS_WARPING_DOFS
        if len(warping) > 0:
            raise NotImplementedError(
                "Element type(s) {} carry warping degrees of freedom, which "
                "are not supported by abqhom.".format(sorted(warping)))

        unsupported = used_types & ELEMENTS_NOT_SUPPORTED
        if len(unsupported) > 0:
            raise NotImplementedError(
                "The element type(s) {} are not supported by abqhom."
                .format(sorted(unsupported)))

        # Initialize Abaqus input file
        path = os.path.abspath(self.model_name + ".inp")
        with open(path, "w") as output:
            output.write("*PART, NAME={}\n".format(self.model_name))

            # Add nodes
            output.write("*NODE\n")
            for tag, coords in zip(node_tags, node_coords):
                output.write("{}, {}, {}, {}\n".format(
                    tag, coords[0], coords[1], coords[2]))

            # Add elements
            for el_type in element_tags:
                type_tags = element_tags[el_type]
                type_conn = conn[el_type]

                output.write(
                    "*ELEMENT, TYPE={}\n".format(el_type_map[el_type]))
                for tag, nodes in zip(type_tags, type_conn):
                    output.write("{}".format(tag))
                    for node in nodes:
                        output.write(", {}".format(node))
                    output.write("\n")

            output.write("*END PART\n")

        return path


class Boundary(UserDict, object):
    """
    Boundary-entity map of an axis-aligned rectangular or cuboidal RVE.

    The class provides a dictionary-like mapping from the boundary-group
    identifiers in :mod:`abqhom.constants` to Gmsh dimension tags. For a
    periodic boundary, primary and secondary entities that correspond under
    a cell translation occupy the same position in their respective lists.

    Parameters
    ----------
    model_name : str
        The name of the underlying Gmsh model.
    initial_data : mapping of str to sequence of (int, int), optional
        Boundary-group names and corresponding Gmsh dimension tags. If None,
        the mapping is generated from the model geometry.
    sampling_density : int, optional
        The number of sampling points per spatial direction that will be
        generated to check for periodicity. The default is 10.
    check_periodic : bool, optional
        If True, check the geometry for periodicity and establish the
        primary--secondary entity correspondence.
    periodic_default : bool, optional
        Periodicity flag used when no automatic check is performed. If True,
        ``initial_data`` must already place corresponding primary and
        secondary dimension tags at the same list index.

    Attributes
    ----------
    periodic : bool
        A flag to mark whether the boundary is periodic.
    model_name : str
        The name of the underlying Gmsh model.
    model_dimension : int
        The geometric dimension of the bounded model.
    sampling_density : int
        The number of sampling points per spatial direction that will be
        generated to check for periodicity.
    data : dict[str, list[tuple[int, int]]]
        Underlying boundary-group mapping.

    """

    def __init__(self, model_name, initial_data=None, sampling_density=10,
                 check_periodic=True, periodic_default=False):

        initial_data = {} if initial_data is None else dict(initial_data)

        if not gmsh.is_initialized():
            gmsh.initialize()
        if model_name not in gmsh.model.list():
            raise ValueError("The Gmsh model {} does not exist.".format(
                model_name))

        super(Boundary, self).__init__(initial_data)
        self._model_name = model_name
        self.sampling_density = sampling_density
        self._periodic = periodic_default

        # Build the boundary automatically
        if not initial_data:
            try:
                self._build_boundary()
            except Exception as error:
                raise RuntimeError(
                    "Automatic boundary generation not possible. Please "
                    "provide the boundary manually. Original error: {}"
                    .format(error))

        # Add missing groups as empty groups
        for key in POINT_GROUPS + EDGE_GROUPS + FACE_GROUPS:
            if key not in self:
                self[key] = []

        # Check if there is any invalid data in the boundary
        if not set(self.keys()) <= set(POINT_GROUPS + EDGE_GROUPS +
                                       FACE_GROUPS):
            raise RuntimeError(
                "The boundary contains data that cannot be associated with "
                "boundary groups.")

        self._validate_group_entities()

        # Rearrange boundary map in case of a periodic geometry
        if check_periodic:
            try:
                self._check_for_periodicity()
            except Exception as error:
                raise RuntimeError(
                    "Automatic periodicity check not possible. Please set "
                    "check_periodic = False and provide information on "
                    "periodicity through construction of initial_data and "
                    "periodic_default. Original error: {}".format(error))

    @property
    def periodic(self):
        """
        A flag to mark whether the boundary is periodic.

        Returns
        -------
        bool
            True if the boundary is periodic, False otherwise.

        """

        return self._periodic

    @property
    def model_name(self):
        """
        The name of the model the boundary belongs to.

        Returns
        -------
        str
            The model name.

        """

        return self._model_name

    @property
    def sampling_density(self):
        """
        Sample points per direction to use for the periodicity check.

        Returns
        -------
        int
            The number of sample points per spatial direction.

        """

        return self._sampling_density

    @sampling_density.setter
    def sampling_density(self, n):
        """
        Sample points per direction to use for the periodicity check.

        Parameters
        ----------
        n : int
            The number of sample points per spatial direction.

        """

        if not isinstance(n, int) or n < 1:
            raise ValueError(
                "The sampling density must be a positive integer.")
        self._sampling_density = int(n)

    @staticmethod
    def _active_axis_mask(bounds_min, bounds_max, tolerance):
        """Return the non-degenerate axes of axis-aligned bounds."""

        spans = np.asarray(bounds_max) - np.asarray(bounds_min)
        return np.invert(
            np.isclose(spans, 0.0, atol=tolerance, rtol=0.0))

    @property
    def model_dimension(self):
        """
        The geometric dimension of the bounded model.

        Returns
        -------
        int
            The number of non-degenerate coordinate axes.

        """

        bounding_box = np.asarray(
            self.get_bounding_box_from_points(),
            dtype=np.float64).reshape(2, 3)
        tolerance = gmsh.option.get_number("Geometry.Tolerance")
        return int(self._active_axis_mask(
            bounding_box[0], bounding_box[1], tolerance).sum())

    def __setitem__(self, key, value):
        if (not isinstance(key, str) or
                key.upper() not in set(POINT_GROUPS +
                                       EDGE_GROUPS +
                                       FACE_GROUPS)):
            raise KeyError("Key must be a group identifier string.")
        super(Boundary, self).__setitem__(key.upper(), value)

    def __getitem__(self, key):
        if (not isinstance(key, str) or
                key.upper() not in set(POINT_GROUPS +
                                       EDGE_GROUPS +
                                       FACE_GROUPS)):
            raise KeyError("Key must be a group identifier string.")
        return super(Boundary, self).__getitem__(key.upper())

    def __iter__(self):
        """Iterate over boundary group identifiers."""

        return iter(self.data)

    def get_bounding_box_from_points(self):
        """
        Return the model bounding box from CAD points or mesh nodes.

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

        # Point coordinates avoid the tolerance expansion applied by
        # gmsh.model.get_bounding_box(), particularly for OCC geometries.

        gmsh.model.set_current(self.model_name)

        # Prefer connected point entities. Isolated construction points do
        # not form part of the RVE support and must not enlarge its bounds.
        points = gmsh.model.get_entities(0)
        connected_points = []
        for point in points:
            upward, _ = gmsh.model.get_adjacencies(*point)
            if len(upward) > 0:
                connected_points.append(point)
        if len(connected_points) > 0:
            points = connected_points

        # Fall back to mesh nodes if no connected CAD points are available,
        # e.g. for a pre-built discrete mesh.
        if len(connected_points) == 0:
            _, node_coords, _ = gmsh.model.mesh.get_nodes()
            if len(node_coords) > 0:
                node_coords = np.array(node_coords).reshape(-1, 3)
                mn, mx = node_coords.min(axis=0), node_coords.max(axis=0)
                return tuple(float(value) for value in list(mn) + list(mx))
            if len(points) == 0:
                return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

        coordinates = np.asarray([
            gmsh.model.get_value(dim, tag, []) for dim, tag in points],
            dtype=np.float64)
        bounds_min = coordinates.min(axis=0)
        bounds_max = coordinates.max(axis=0)
        return tuple(float(value) for value in
                     np.concatenate((bounds_min, bounds_max)))

    def _validate_group_entities(self):
        """Validate the entity references in the boundary mapping."""

        gmsh.model.set_current(self.model_name)
        available_entities = set(gmsh.model.get_entities())
        assigned_entities = {}

        for group_name, dimtags in self.items():
            for dimtag in dimtags:
                if (not isinstance(dimtag, (tuple, list))
                        or len(dimtag) != 2):
                    raise RuntimeError(
                        "Boundary group '{}' contains an invalid Gmsh "
                        "dimension tag: {}."
                        .format(group_name, dimtag))
                entity = (int(dimtag[0]), int(dimtag[1]))
                if entity not in available_entities:
                    raise RuntimeError(
                        "Boundary group '{}' references a Gmsh entity that "
                        "does not exist: {}."
                        .format(group_name, entity))
                if entity in assigned_entities:
                    raise RuntimeError(
                        "Gmsh entity {} is assigned to both '{}' and '{}'."
                        .format(entity, assigned_entities[entity],
                                group_name))
                assigned_entities[entity] = group_name

    def _build_boundary(self):
        """Build the boundary-entity map from the model geometry."""

        # Activate the model
        gmsh.model.set_current(self.model_name)

        # Get the exact cell bounds from its points. Entity bounding boxes
        # can contain a small kernel-dependent tolerance expansion.
        bounding_box = np.asarray(
            self.get_bounding_box_from_points(),
            dtype=np.float64).reshape(2, 3)
        bounds_min, bounds_max = bounding_box
        geometry_tol = gmsh.option.get_number("Geometry.Tolerance")
        group_to_coords = self.group_coordinate_identifiers(
            bounds_min, bounds_max, tolerance=geometry_tol)
        coordinate_scale = max(
            1.0, float(np.max(np.abs(bounding_box))))
        roundoff = 10.0*np.finfo(np.float64).eps*coordinate_scale
        # OpenCASCADE expands entity boxes by 1e-7 independently of the
        # configurable Gmsh geometry tolerance.
        box_tol = max(geometry_tol, 1.0e-7) + roundoff
        model_dimension = self.model_dimension

        # Initialize boundary map
        for group_name in group_to_coords:
            self[group_name] = []

        # Classify an entity from the coordinates that remain constant over
        # its complete bounding box. This avoids parametrizing irrelevant
        # internal or trimmed entities, which can be unreliable for
        # degenerate OpenCASCADE curves.
        entities = [entity for entity in gmsh.model.get_entities()
                    if entity[0] < 3]
        for dim, tag in entities:
            if dim == 0:
                upward, downward = gmsh.model.get_adjacencies(dim, tag)
                if len(upward) == len(downward) == 0:
                    continue

            entity_box = np.asarray(
                gmsh.model.get_bounding_box(dim, tag), dtype=np.float64)
            entity_min, entity_max = entity_box[:3], entity_box[3:]
            below_bounds = np.any(entity_min < bounds_min - box_tol)
            above_bounds = np.any(entity_max > bounds_max + box_tol)
            if below_bounds or above_bounds:
                raise RuntimeError(
                    "Entity ({}, {}) is outside of the computed "
                    "bounding box.".format(dim, tag))

            # Assign the most specific matching group. Lower-dimensional
            # structural entities may lie in a cell edge or face, so the
            # entity dimension is not required to equal the group dimension.
            for group_name, identifiers in group_to_coords.items():
                group_dimension = int(np.isnan(identifiers).sum())
                if group_dimension >= model_dimension:
                    continue

                mask = np.invert(np.isnan(identifiers))
                matches_min = np.allclose(
                    entity_min[mask], identifiers[mask],
                    atol=box_tol, rtol=0.0)
                matches_max = np.allclose(
                    entity_max[mask], identifiers[mask],
                    atol=box_tol, rtol=0.0)
                if matches_min and matches_max:
                    self[group_name].append((dim, tag))
                    break

        # Cleanup of boundary groups
        for group_name, group_dimtags in self.items():
            to_remove = set()
            group_dimtag_set = set(group_dimtags)
            for dimtag in group_dimtags:
                _, downward = gmsh.model.get_adjacencies(*dimtag)

                # Remove entities represented through the boundary of a
                # higher-dimensional entity in the same group.
                downward = {(dimtag[0] - 1, tag) for tag in downward}
                to_remove.update(downward & group_dimtag_set)

            self[group_name] = [dimtag for dimtag in group_dimtags
                                if dimtag not in to_remove]

    def get_group_coordinate_identifiers(self):
        """
        Return boundary-group identifiers for the current model bounds.

        Returns
        -------
        group_to_coords : OrderedDict[str, ndarray]
            The group identifier map. Keys are the group names and
            values are the identifier coordinates. Coordinates that are free
            along an edge or face are marked with ``numpy.nan``.

        """

        bounding_box = np.asarray(
            self.get_bounding_box_from_points(),
            dtype=np.float64).reshape(2, 3)
        return self.group_coordinate_identifiers(
            bounding_box[0], bounding_box[1],
            tolerance=gmsh.option.get_number("Geometry.Tolerance"))

    @staticmethod
    def group_coordinate_identifiers(bounds_min, bounds_max, tolerance=0.0):
        """
        Return boundary-group identifiers for supplied axis-aligned bounds.

        Two-dimensional bounds use the canonical groups
        ``POINT_1``--``POINT_4`` and ``EDGE_1``--``EDGE_4``.

        Parameters
        ----------
        bounds_min : array_like
            Minimum coordinates ``(x_min, y_min, z_min)``.
        bounds_max : array_like
            Maximum coordinates ``(x_max, y_max, z_max)``.
        tolerance : float, optional
            Absolute tolerance used to identify degenerate coordinate axes.
            The default is zero.

        Returns
        -------
        group_to_coords : OrderedDict[str, ndarray]
            The group identifier map. Keys are the group names and values
            are the identifier coordinates. Coordinates that are free along
            an edge or face are marked with ``numpy.nan``.

        """

        try:
            bounds_min = np.asarray(bounds_min, dtype=np.float64)
            bounds_max = np.asarray(bounds_max, dtype=np.float64)
            tolerance = float(tolerance)
        except (TypeError, ValueError):
            raise ValueError(
                "Boundary bounds and tolerance must be numeric.")

        if bounds_min.shape != (3,) or bounds_max.shape != (3,):
            raise ValueError("Boundary bounds must have shape (3,).")
        minimum_is_finite = np.all(np.isfinite(bounds_min))
        maximum_is_finite = np.all(np.isfinite(bounds_max))
        if not minimum_is_finite or not maximum_is_finite:
            raise ValueError("Boundary bounds must be finite.")
        if np.any(bounds_min > bounds_max):
            raise ValueError(
                "Each minimum boundary coordinate must not exceed its "
                "maximum.")
        if not np.isfinite(tolerance) or tolerance < 0.0:
            raise ValueError("The boundary tolerance must be nonnegative.")

        active_axes = np.flatnonzero(Boundary._active_axis_mask(
            bounds_min, bounds_max, tolerance))
        dimension = len(active_axes)
        group_identifiers = np.full((26, 3), np.nan, dtype=np.float64)

        # Signature entries -1, 0, and 1 denote a coordinate at the
        # minimum, free along the entity, and at the maximum, respectively.
        # Their order follows POINT_GROUPS, EDGE_GROUPS, and FACE_GROUPS.
        signatures_3d = [
            (-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1),
            (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1),
            (0, -1, -1), (1, 0, -1), (0, 1, -1), (-1, 0, -1),
            (0, -1, 1), (1, 0, 1), (0, 1, 1), (-1, 0, 1),
            (-1, -1, 0), (1, -1, 0), (1, 1, 0), (-1, 1, 0),
            (0, 0, -1), (0, 0, 1), (0, -1, 0), (0, 1, 0),
            (-1, 0, 0), (1, 0, 0)]
        signatures_2d = [(-1, -1), (1, -1), (1, 1), (-1, 1),
                         (0, -1), (1, 0), (0, 1), (-1, 0)]

        if dimension == 3:
            rows = range(26)
            signatures = signatures_3d
        elif dimension == 2:
            rows = list(range(4)) + list(range(8, 12))
            signatures = signatures_2d
        elif dimension == 1:
            rows = [0, 1]
            signatures = [(-1,), (1,)]
        else:
            rows = [0]
            signatures = [()]

        for row, signature in zip(rows, signatures):
            identifiers = np.array(bounds_min, dtype=np.float64)
            for axis, location in zip(active_axes, signature):
                if location < 0:
                    identifiers[axis] = bounds_min[axis]
                elif location > 0:
                    identifiers[axis] = bounds_max[axis]
                else:
                    identifiers[axis] = np.nan
            group_identifiers[row] = identifiers

        # Assign the identifiers to the group names
        group_to_coords = OrderedDict()
        for i, j in enumerate(POINT_GROUPS + EDGE_GROUPS + FACE_GROUPS):
            group_to_coords[j] = group_identifiers[i]

        return group_to_coords

    @staticmethod
    def _generate_sample_points_on_entity(dim, tag, n, tol):
        """Sample an entity in each of its parametric directions."""

        (bounds_min,
         bounds_max) = gmsh.model.get_parametrization_bounds(dim, tag)
        bounds_min = np.asarray(bounds_min, dtype=np.float64)
        bounds_max = np.asarray(bounds_max, dtype=np.float64)
        sample_points = [np.linspace(i, j, n)
                         for i, j in zip(bounds_min + tol, bounds_max - tol)]
        sample_points = np.meshgrid(*sample_points)

        # Get the physical coordinates from the parametric space
        assert len(sample_points) in {0, 1, 2}
        sample_coords = []  # Points
        if len(sample_points) == 1:  # Lines
            sample_coords = sample_points[0]
        elif len(sample_points) == 2:  # Faces
            u, v = sample_points
            sample_coords = np.vstack((u.flatten(), v.flatten())).T.flatten()
        coords = gmsh.model.get_value(dim, tag, sample_coords)
        coords = np.array(coords, dtype=np.float64).reshape(-1, 3)

        # Filter sample points of trimmed entities to make sure to only
        # obtain points inside the trimmed domain
        if dim > 0:
            inside = np.array([gmsh.model.is_inside(dim, tag, c)
                               for c in coords])
            coords = coords[inside.nonzero()]

        return coords

    def _check_for_periodicity(self):
        """Check and pair translationally periodic boundary entities."""

        group_to_coords = self.get_group_coordinate_identifiers()

        def sort_sample_coordinates(coords, primary):
            """Sort sampled coordinates along free entity directions."""

            order = np.where(np.isnan(group_to_coords[primary]))[0][::-1]
            order = order if len(order) > 0 else [0]
            rounded_coords = np.round(coords, decimals=10)
            coords = coords[np.lexsort([rounded_coords[:, i] for i in order])]
            return coords

        boundary_groups = {0: POINT_GROUPS[0:1],
                           1: POINT_GROUPS[0:2],
                           2: POINT_GROUPS[0:4] + EDGE_GROUPS[0:4],
                           3: POINT_GROUPS + EDGE_GROUPS + FACE_GROUPS}

        # Make a deepcopy of the original boundary map
        original_boundary_map = {i: list(j) for i, j in self.data.items()}

        self._periodic = True
        tol = gmsh.option.get_number("Geometry.Tolerance")
        model_dimension = self.model_dimension
        for secondary, primary in PERIODIC_PAIRS.items():
            secondary_indices = []
            used_secondary_indices = set()
            primary_entities = self[primary]
            secondary_entities = self[secondary]

            # Skip irrelevant secondary groups
            if secondary not in boundary_groups[model_dimension]:
                continue

            # A periodic entity correspondence must be bijective. Check the
            # original group cardinalities before replacing the secondary
            # list with its matched order.
            if len(primary_entities) != len(secondary_entities):
                self._periodic = False
                self.data = original_boundary_map
                return

            for primary_dimtag in primary_entities:
                # Generate sample points on the primary entity
                primary_coords = self._generate_sample_points_on_entity(
                    primary_dimtag[0], primary_dimtag[1],
                    self.sampling_density,
                    tol)
                if len(primary_coords) == 0:
                    continue
                primary_coords = sort_sample_coordinates(primary_coords,
                                                         primary)

                # Check periodicity for all potential secondary entities
                for secondary_index, secondary_dimtag in enumerate(
                        secondary_entities):
                    if secondary_index in used_secondary_indices:
                        continue
                    if primary_dimtag[0] != secondary_dimtag[0]:
                        continue

                    # Generate sample points on the potential secondary
                    # entity
                    secondary_coords = self._generate_sample_points_on_entity(
                        secondary_dimtag[0], secondary_dimtag[1],
                        self.sampling_density, tol)
                    if len(secondary_coords) == 0:
                        continue
                    secondary_coords = sort_sample_coordinates(
                        secondary_coords, primary)

                    # If the number of points generated is not the same
                    # the entities cannot be periodically connected
                    if len(primary_coords) != len(secondary_coords):
                        continue

                    # Check if the primary and secondary entities are
                    # periodically connected
                    delta_true = np.abs(secondary_coords - primary_coords)
                    delta_expected = np.nan_to_num(np.abs(
                        group_to_coords[secondary]
                        - group_to_coords[primary]))

                    if np.allclose(delta_true, delta_expected, atol=tol,
                                   rtol=0.0):
                        secondary_indices.append(secondary_index)
                        used_secondary_indices.add(secondary_index)
                        break

            # Update the order of the secondary entities and check
            # periodicity
            self[secondary] = [secondary_entities[i]
                               for i in secondary_indices]
            if len(self[secondary]) != len(primary_entities):
                self._periodic = False
                self.data = original_boundary_map
                return
