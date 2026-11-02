# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Module containing the AbaqusRVE class that provides useful methods for
computational homogenization with Abaqus.

"""

import os
from numbers import Integral

import numpy as np

from abqhom.constants import PERIODIC_PAIRS, ELEMENTS_ROTATIONAL_DOFS

from abaqus import Mdb, mdb, session
from regionToolset import Region
from mesh import MeshNodeArray, MeshElementArray
from abaqusConstants import (THREE_D, TWO_D_PLANAR, DEFORMABLE_BODY, ON, OFF,
                             UNSET, CARTESIAN, N1_COSINES, DURING_ANALYSIS,
                             BEFORE_ANALYSIS, TOTAL, TIME, COMPLETED)
from paraqus.abaqus import OdbReader
from paraqus.writers import BinaryWriter, CollectionWriter


class AbaqusRVE(object):
    """
    Abaqus representation of a periodic structural RVE.

    The class imports a :class:`GmshRVE`, transfers its physical and boundary
    sets, applies macroscopic loading and periodic constraints, and recovers
    effective quantities from the Abaqus output database. By default,
    redundant element copies on secondary periodic boundary entities are
    removed while all nodes are retained so that periodic constraints can
    subsequently be applied. The resulting discretization is an
    ownership-consistent reduced model.

    Constructing an instance initializes a new Abaqus model database through
    ``Mdb()`` and therefore replaces the active in-memory model database.

    Parameters
    ----------
    gmsh_RVE : GmshRVE
        Gmsh representation containing the RVE mesh and boundary sets.
    custom_types : mapping of int to str, optional
        Overrides for the default Gmsh-to-Abaqus element-type map.
        Element types with warping degrees of freedom are not supported.
    remove_secondary : bool, optional
        If True, redundant periodic element copies are removed. Setting
        this to False retains the complete periodic representation for
        inspection or visualization, but the resulting unweighted problem
        is not ownership consistent. The default is True.

    Attributes
    ----------
    model_name : str
        The name of the model.
    periodic_mesh : bool
        A flag indicating whether opposite boundary node sets form
        one-to-one periodic pairs.
    model_dimension : SymbolicConstant
        The model dimension. The value can be ``THREE_D`` or
        ``TWO_D_PLANAR``.
    model : Model
        The internal Abaqus representation of the model.
    part : Part
        The internal Abaqus representation of the RVE part.
    instance : PartInstance
        The internal representation of the RVE part in the assembly.
    odb : Odb or None
        The Abaqus output database containing the analysis results, or None
        when no database is open.
    job_name : str
        Name of the most recently submitted job, or the default name before
        the first analysis.

    """

    def __init__(self, gmsh_RVE, custom_types=None, remove_secondary=True):
        if gmsh_RVE.model_dimension not in {2, 3}:
            raise RuntimeError("Only 2D and 3D models are supported.")
        dim = THREE_D if gmsh_RVE.model_dimension == 3 else TWO_D_PLANAR
        custom_types = {} if custom_types is None else dict(custom_types)

        # Write Abaqus input file
        input_file = gmsh_RVE.write_abaqus_input(custom_types=custom_types)

        # Load the RVE into Abaqus
        Mdb()
        default_model = mdb.models.keys()[0]
        mdb.ModelFromInputFile(name=gmsh_RVE.model_name,
                               inputFileName=input_file)
        del mdb.models[default_model]

        self._model_name = gmsh_RVE.model_name
        self._periodic_mesh = gmsh_RVE.periodic_mesh
        self._model_dimension = dim
        self._bc_initialized = False
        self._has_periodic_bc = False
        self._model = mdb.models[self._model_name]
        self._part = self._model.parts[self._model_name]
        self._instance = None
        self._odb = None
        self._job_name = None
        self._job_count = 0
        self._section_radius = {}

        # Store the geometric bounding box of the RVE
        self._bounding_box = tuple(
            float(v) for v in gmsh_RVE.get_bounding_box_from_points())

        # Set the correct modeling space
        self._part.setValues(space=dim)

        # Add node sets
        if "ALL_NODES" not in gmsh_RVE.node_sets:
            gmsh_RVE.build_node_sets()
        self._add_node_sets(gmsh_RVE.node_sets)

        # Add element sets
        if "ALL_ELEMENTS" not in gmsh_RVE.element_sets:
            gmsh_RVE.build_element_sets()
        self._add_element_sets(gmsh_RVE.element_sets)

        # Record the element-supported nodal DOFs before periodic element
        # copies are removed. This retains the original structural trace for
        # the subsequent construction of periodic rotational constraints.
        self._unreduced_element_types = self._collect_element_types()
        self._unreduced_node_element_types = (
            self._collect_node_element_types())

        # Remove redundant periodic copies while retaining their nodes
        if (remove_secondary and self.periodic_mesh
                and "ELSET_SECONDARY_ELEMENTS" in self.part.sets):
            self.part.deleteElement(
                elements=self.part.sets["ELSET_SECONDARY_ELEMENTS"],
                deleteUnreferencedNodes=OFF)

        # Add part instance to the assembly
        assembly = self.model.rootAssembly
        assembly.DatumCsysByDefault(CARTESIAN)
        self._instance = assembly.Instance(name=self._model_name,
                                           part=self._part,
                                           dependent=ON)

    def __str__(self):
        dim = 3 if self.model_dimension == THREE_D else 2
        descr = "AbaqusRVE '{}'\n".format(self.model_name)
        descr += "\tmodel dimension       : {}\n".format(dim)
        descr += "\t# steps               : {}\n".format(
            len(self.model.steps))
        descr += "\t# constraints         : {}\n".format(
            len(self.model.constraints))
        descr += "\t# boundary conditions : {}\n".format(
            len(self.model.boundaryConditions))
        descr += "\t# load conditions     : {}".format(len(self.model.loads))
        return descr

    def __repr__(self):
        return self.__str__()

    @property
    def model_name(self):
        """
        The internal name of the Abaqus model.

        Returns
        -------
        str
            The model name.

        """

        return self._model_name

    @property
    def periodic_mesh(self):
        """
        Whether opposite boundary node sets form periodic pairs.

        Returns
        -------
        bool
            True if the boundary nodes can be paired one to one, False
            otherwise.

        """

        return self._periodic_mesh

    @property
    def model_dimension(self):
        """
        The geometric dimension of the Abaqus model.

        Returns
        -------
        SymbolicConstant
            The dimension, i.e. THREE_D or TWO_D_PLANAR.

        """

        return self._model_dimension

    @property
    def model(self):
        """
        The internal Abaqus representation of the model.

        Returns
        -------
        Model
            The Abaqus model.

        """

        return self._model

    @property
    def part(self):
        """
        The internal Abaqus representation of the RVE part.

        Returns
        -------
        Part
            The RVE Abaqus part.

        """

        return self._part

    @property
    def instance(self):
        """
        The internal Abaqus representation of the RVE part instance.

        Returns
        -------
        PartInstance
            The RVE Abaqus part instance.

        """

        return self._instance

    @property
    def odb(self):
        """
        The Abaqus output database containing the analysis results.

        Returns
        -------
        Odb or None
            The open output database, or None.

        """
        return self._odb

    @property
    def job_name(self):
        """
        Name of the most recently submitted Abaqus job.

        Before the first analysis, the default name ``JOB_<model>`` is
        returned. Repeated analyses receive a numeric suffix to avoid reuse
        of solver resources associated with an earlier job.

        Returns
        -------
        str
            The current Abaqus job name.

        """

        name = getattr(self, "_job_name", None)
        return name if name is not None else "JOB_" + self.model_name

    def _add_node_sets(self, node_sets):
        """Add the supplied node sets to the Abaqus model."""

        for set_name, set_items in node_sets.items():
            set_name = "NSET_" + str(set_name)
            set_items = [int(i) for i in set_items]

            # Sets must be unsorted to not break ordering for periodic BCs
            if len(set_items) > 0:
                self.part.SetFromNodeLabels(name=set_name,
                                            nodeLabels=set_items,
                                            unsorted=True)

    def _add_element_sets(self, element_sets):
        """Add the supplied element sets to the Abaqus model."""

        for set_name, set_items in element_sets.items():
            set_name = "ELSET_" + str(set_name)
            set_items = [int(i) for i in set_items]

            if len(set_items) > 0:
                self.part.SetFromElementLabels(name=set_name,
                                               elementLabels=set_items)

    def get_bounding_box(self):
        """
        Return the geometric bounding box of the RVE.

        The bounds are stored at construction from the Gmsh geometry or,
        for a pre-built discrete mesh, from its nodal coordinates.

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

        x_min, y_min, z_min, x_max, y_max, z_max = self._bounding_box
        return x_min, y_min, z_min, x_max, y_max, z_max

    def get_volume(self):
        r"""
        Return the reference-cell measure :math:`V=|\Omega_0|`.

        Returns
        -------
        float
            Reference-cell volume in 3D or area per unit out-of-plane
            thickness in 2D.

        """

        bounding_box = np.array(self.get_bounding_box()).reshape(2, -1)
        delta = bounding_box[1] - bounding_box[0]
        if self.model_dimension == TWO_D_PLANAR:
            delta[-1] = 1.0
        return float(np.prod(delta))

    def _collect_element_types(self):
        """Return the names of all element types in the part mesh."""

        all_elements = self.part.sets["ELSET_ALL_ELEMENTS"].elements
        element_types = {str(element.type) for element in all_elements}
        return element_types

    def _collect_node_element_types(self):
        """Map each node label to the types of its attached elements."""

        node_element_types = {}
        for element in self.part.elements:
            element_type = str(element.type)
            for index in element.connectivity:
                label = self.part.nodes[index].label
                node_element_types.setdefault(label, set()).add(element_type)
        return node_element_types

    def _embedded_region_nodes(self):
        """Return labels of nodes tied by embedded-region constraints."""

        labels = set()
        for name in self.model.constraints.keys():
            constraint = self.model.constraints[name]
            if constraint.__class__.__name__ != "EmbeddedRegion":
                continue

            set_name = constraint.embeddedRegion[0]
            if set_name in self.part.sets:
                labels.update(n.label for n in self.part.sets[set_name].nodes)

        return labels

    def _add_dirichlet_from_displacement_history(self, time_history, u_history,
                                                 region, step_name, bc_name,
                                                 amplitude_name):
        """Add history-dependent displacement boundary conditions."""

        for i in range(3):
            # Ignore the out-of-plane component for 2D models
            if i == 2 and self.model_dimension != THREE_D:
                break

            # Add amplitude for the load history
            amplitude_name_i = amplitude_name + "_{}".format(i + 1)
            data = np.vstack((time_history, u_history[:, i])).T
            data = tuple([(t, v) for t, v in data])
            self.model.TabularAmplitude(name=amplitude_name_i, timeSpan=TOTAL,
                                        data=data)

            # Add the boundary condition
            # The actual displacement value will be controlled by the
            # amplitude, so just put a value of 1.0 here
            bc_name_i = bc_name + "_{}".format(i + 1)
            u = [UNSET, UNSET, UNSET]
            u[i] = 1.0
            self.model.DisplacementBC(name=bc_name_i,
                                      createStepName=step_name,
                                      region=region,
                                      amplitude=amplitude_name_i,
                                      u1=u[0], u2=u[1], u3=u[2])

    def _clear_boundary_conditions(self):
        """Remove boundary conditions and constraints created by this class."""

        model = self.model
        assembly = model.rootAssembly

        for name in list(model.boundaryConditions.keys()):
            if name.startswith("DIRICHLET"):
                del model.boundaryConditions[name]
        for name in list(model.amplitudes.keys()):
            if name.startswith("DIRICHLET"):
                del model.amplitudes[name]
        for name in list(model.constraints.keys()):
            if name.startswith("PBC_"):
                del model.constraints[name]
        for name in list(assembly.instances.keys()):
            if name.startswith("REFERENCE_POINT_"):
                del assembly.instances[name]
        for name in list(model.parts.keys()):
            if name.startswith("REFERENCE_POINT_"):
                del model.parts[name]

        self._bc_initialized = False
        self._has_periodic_bc = False

    def apply_affine_displacement_BC(self, initialization_step,
                                     displacement_gradient_history):
        r"""
        Apply affine displacement boundary conditions.

        Every boundary node is prescribed according to

        .. math::
            \boldsymbol{u}(\boldsymbol{X},\tau)
            =\bar{\boldsymbol{H}}(\tau)\boldsymbol{X},

        where :math:`\bar{\boldsymbol{H}}=\bar{\boldsymbol{F}}-\boldsymbol{I}`.
        For geometrically linear strain loading, the skew part is set to zero
        and :math:`\bar{\boldsymbol{H}}=\bar{\boldsymbol{\varepsilon}}`.
        The history uses total analysis time, not step time, because the
        boundary conditions propagate to subsequent steps.

        **Elements with rotational degrees of freedom, such as beams and
        shells, are not supported by the affine displacement boundary
        conditions.**

        Parameters
        ----------
        initialization_step : str
            The name of the step in which the boundary conditions are
            initialized.
        displacement_gradient_history : LoadHistory
            History of the macroscopic displacement gradient
            :math:`\bar{\boldsymbol{H}}(\tau)`.

        Returns
        -------
        None

        """

        if self._bc_initialized:
            self._clear_boundary_conditions()

        # Rotational DOFs (beams and shells) are not supported: the
        # consistent boundary rotation depends on the macroscopic
        # rotation, which is only defined once a macroscopic
        # rotation/curvature model is fixed. Use periodic BCs instead.
        if len(self._collect_element_types() & ELEMENTS_ROTATIONAL_DOFS) > 0:
            raise NotImplementedError(
                "Affine displacement boundary conditions are not supported "
                "for elements with rotational degrees of freedom "
                "(beams and shells). Please use periodic "
                "boundary conditions instead.")

        time_history = displacement_gradient_history.time_data
        displacement_gradients = displacement_gradient_history.load_data

        # Apply Dirichlet boundary conditions on all boundary nodes
        boundary_nodes = self.instance.sets["NSET_BOUNDARY_NODES"]
        for node in boundary_nodes.nodes:
            region = Region(nodes=MeshNodeArray(nodes=[node]))
            x = np.array(node.coordinates)
            u_history = np.dot(displacement_gradients, x)

            amplitude_name = "DIRICHLET_{}".format(node.label)
            bc_name = "DIRICHLET_{}".format(node.label)
            self._add_dirichlet_from_displacement_history(
                time_history, u_history, region, initialization_step, bc_name,
                amplitude_name)

        self._bc_initialized = True

    def apply_periodic_BC(self, initialization_step,
                          displacement_gradient_history, fixed_node=None):
        r"""
        Apply periodic translational and rotational constraints.

        For every primary--secondary node pair :math:`(p,s)` separated by
        the cell translation :math:`\Delta\boldsymbol{X}^{g}`, the method
        enforces

        .. math::
            \boldsymbol{u}_{s}(\tau)-\boldsymbol{u}_{p}(\tau)
            =\bar{\boldsymbol{H}}(\tau)\Delta\boldsymbol{X}^{g}.

        Here :math:`\bar{\boldsymbol{H}}=\bar{\boldsymbol{F}}-\boldsymbol{I}`,
        with :math:`\bar{\boldsymbol{H}}=\bar{\boldsymbol{\varepsilon}}` for
        geometrically linear strain loading. Nodes carrying rotational
        degrees of freedom additionally satisfy

        .. math::
            \boldsymbol{\phi}_{s}(\tau)-\boldsymbol{\phi}_{p}(\tau)
            =\boldsymbol{0},

        which realizes periodic microrotation with zero imposed macroscopic
        curvature. The history uses total analysis time, not step time.

        Nodes whose translational degrees of freedom are eliminated by an
        embedded-region constraint (e.g. fibers embedded in a matrix) are
        excluded from the translational equations, since Abaqus cannot form
        an equation on an eliminated DOF. Their translational periodicity
        is inherited from the periodic host. Any rotational DOFs are still
        coupled as above.

        Parameters
        ----------
        initialization_step : str
            The name of the step in which the boundary conditions are
            initialized.
        displacement_gradient_history : LoadHistory
            History of the macroscopic displacement gradient
            :math:`\bar{\boldsymbol{H}}(\tau)`.
        fixed_node : str, optional
            The name of a node set whose first admissible node is used
            as the anchor. If None, an anchor is selected automatically.
            The default is None.

        Returns
        -------
        None

        """

        if self._bc_initialized:
            self._clear_boundary_conditions()

        time_history = displacement_gradient_history.time_data
        displacement_gradients = displacement_gradient_history.load_data

        def create_reference_point(primary, secondary):
            """Create the control point for a periodic entity pair."""

            coords = [0.0, 0.0, 0.0]
            point_name = "REFERENCE_POINT_{}_{}".format(primary, secondary)
            part = self.model.Part(dimensionality=self.model_dimension,
                                   name=point_name, type=DEFORMABLE_BODY)
            point = self.model.parts[point_name].ReferencePoint(point=coords)
            self.model.parts[point_name].Set(
                name="REFERENCE_POINT",
                referencePoints=[self.model.parts[point_name]
                                 .referencePoints[point.id]])
            self.model.rootAssembly.Instance(name=point_name, part=part,
                                             dependent=ON)

            # Compute the +/- coordinates of the respective groups
            x_minus = RVE_sets["NSET_" + primary].nodes[0].coordinates
            x_plus = RVE_sets["NSET_" + secondary].nodes[0].coordinates

            # Apply the boundary conditions
            region = self.model.rootAssembly.instances[point_name].sets[
                "REFERENCE_POINT"]
            u_history = np.dot(displacement_gradients,
                               np.array(x_plus) - np.array(x_minus))

            amplitude_name = "DIRICHLET_{}".format(point_name)
            bc_name = "DIRICHLET_{}".format(point_name)
            self._add_dirichlet_from_displacement_history(
                time_history, u_history, region, initialization_step, bc_name,
                amplitude_name)

        if not self.periodic_mesh:
            raise RuntimeError("The mesh is not suitable for periodic BCs.")

        # Define some shortcuts
        RVE_sets = self.instance.sets
        dim = 3 if self.model_dimension == THREE_D else 2

        # Suppress rigid-body translation
        region, anchor_node = self.select_rigid_body_anchor(fixed_node)
        x = np.array(anchor_node.coordinates)
        u_history = np.dot(displacement_gradients, x)
        amplitude_name = "DIRICHLET_ANCHOR_{}".format(anchor_node.label)
        bc_name = "DIRICHLET_ANCHOR_{}".format(anchor_node.label)
        self._add_dirichlet_from_displacement_history(
            time_history, u_history, region, initialization_step, bc_name,
            amplitude_name)

        # Nodes tied by an embedded-region constraint
        embedded = self._embedded_region_nodes()

        # Add periodic coupling to edges and faces
        for secondary, primary in PERIODIC_PAIRS.items():
            secondary_set = "NSET_" + secondary
            primary_set = "NSET_" + primary
            if primary_set not in RVE_sets or secondary_set not in RVE_sets:
                continue

            # Restrict the equation to node pairs whose translations are not
            # eliminated by an embedded constraint
            if embedded:
                secondary_labels, primary_labels = [], []
                secondary_nodes = self.part.sets[secondary_set].nodes
                primary_nodes = self.part.sets[primary_set].nodes
                if len(secondary_nodes) != len(primary_nodes):
                    raise RuntimeError(
                        "Nodes of periodic sets '{}' and '{}' do not form "
                        "one-to-one pairs.".format(secondary_set,
                                                   primary_set))

                # Retain independent pairs and reject inconsistent pairs
                for secondary_node, primary_node in zip(
                        secondary_nodes, primary_nodes):

                    secondary_embedded = secondary_node.label in embedded
                    primary_embedded = primary_node.label in embedded
                    if secondary_embedded != primary_embedded:
                        raise RuntimeError(
                            "Periodic node pair ({}, {}) is embedded on only "
                            "one side of the RVE.".format(
                                primary_node.label, secondary_node.label))

                    if not secondary_embedded:
                        secondary_labels.append(secondary_node.label)
                        primary_labels.append(primary_node.label)

                if not secondary_labels:
                    continue  # all pairs are embedded and handled by the host

                # Create subset of nodes containing only the nodes not
                # handled by the embedded region host
                if (len(secondary_labels)
                        < len(self.part.sets[secondary_set].nodes)):
                    secondary_set = "NSET_{}_INDEPENDENT".format(secondary)
                    primary_set = "NSET_{}_INDEPENDENT".format(primary)

                    if secondary_set not in self.part.sets:
                        self.part.SetFromNodeLabels(
                            name=secondary_set,
                            nodeLabels=secondary_labels,
                            unsorted=True)
                    if primary_set not in self.part.sets:
                        self.part.SetFromNodeLabels(
                            name=primary_set,
                            nodeLabels=primary_labels,
                            unsorted=True)

            # Create a reference point for the right-hand side
            create_reference_point(primary, secondary)

            base_name = self.model_name + "."
            for dof in range(1, dim + 1):
                ref_set_name = ("REFERENCE_POINT_{}_{}.REFERENCE_POINT"
                                .format(primary, secondary))
                terms = [[1.0, base_name + secondary_set, dof],
                         [-1.0, base_name + primary_set, dof],
                         [-1.0, ref_set_name, dof]]

                self.model.Equation(name="PBC_{}_{}_{}".format(
                    primary, secondary, dof), terms=terms)

        # Apply constraints also to rotations
        self._periodic_boundary_rotations()

        self._bc_initialized = True
        self._has_periodic_bc = True

    def _periodic_boundary_rotations(self):
        """Couple nodal rotations periodically using unreduced connectivity."""

        if not (self._unreduced_element_types & ELEMENTS_ROTATIONAL_DOFS):
            return

        rotational_dofs = ([4, 5, 6] if self.model_dimension == THREE_D
                           else [6])

        node_element_types = self._unreduced_node_element_types

        # Use the unreduced structural trace so that a node orphaned by the
        # ownership reduction remains periodically tied to its partner. Pure
        # continuum pairs remain excluded and acquire no artificial rotations.
        for secondary, primary in PERIODIC_PAIRS.items():
            secondary_set = "NSET_" + secondary
            primary_set = "NSET_" + primary
            if (primary_set not in self.instance.sets
                    or secondary_set not in self.instance.sets):
                continue

            for dof in rotational_dofs:
                secondary_nodes = self.part.sets[secondary_set].nodes
                primary_nodes = self.part.sets[primary_set].nodes
                if len(secondary_nodes) != len(primary_nodes):
                    raise RuntimeError(
                        "Nodes of periodic sets '{}' and '{}' do not form "
                        "one-to-one pairs.".format(secondary_set,
                                                   primary_set))

                # Retain pairs for which both nodes possessed the respective
                # DOF before redundant periodic element copies were removed.
                paired_labels = []
                for secondary_node, primary_node in zip(secondary_nodes,
                                                        primary_nodes):
                    secondary_types = node_element_types.get(
                        secondary_node.label, set())
                    primary_types = node_element_types.get(
                        primary_node.label, set())
                    if ((secondary_types & ELEMENTS_ROTATIONAL_DOFS)
                            and (primary_types & ELEMENTS_ROTATIONAL_DOFS)):
                        paired_labels.append(
                            (secondary_node.label, primary_node.label))

                if not paired_labels:
                    continue
                secondary_labels, primary_labels = zip(*paired_labels)

                secondary_subset = "{}_DOF_{}".format(secondary_set, dof)
                primary_subset = "{}_DOF_{}".format(primary_set, dof)
                if secondary_subset not in self.part.sets:
                    self.part.SetFromNodeLabels(
                        name=secondary_subset, nodeLabels=secondary_labels,
                        unsorted=True)
                if primary_subset not in self.part.sets:
                    self.part.SetFromNodeLabels(
                        name=primary_subset, nodeLabels=primary_labels,
                        unsorted=True)

                base_name = self.model_name + "."
                terms = [[1.0, base_name + secondary_subset, dof],
                         [-1.0, base_name + primary_subset, dof]]

                self.model.Equation(name="PBC_{}_{}_{}".format(
                    primary, secondary, dof), terms=terms)

    def select_rigid_body_anchor(self, fixed_node=None):
        """
        Select an independent node for suppressing rigid-body translation.

        The selected node is either an interior node or a node on a primary
        periodic boundary entity. Nodes belonging to an embedded region and
        nodes on secondary periodic entities are excluded.

        Parameters
        ----------
        fixed_node : str, optional
            Node-set name from which the anchor shall be selected. The
            ``NSET_`` prefix is optional. If omitted, a suitable node is
            selected automatically.

        Returns
        -------
        region
            Abaqus region containing the selected node.
        MeshNode
            Selected Abaqus mesh node.

        """

        RVE_sets = self.instance.sets
        embedded_labels = self._embedded_region_nodes()

        # Collect the labels of all periodic secondary nodes
        secondary_labels = set()
        for secondary in PERIODIC_PAIRS.keys():
            set_name = "NSET_" + secondary
            if set_name in RVE_sets:
                secondary_labels.update(
                    n.label for n in RVE_sets[set_name].nodes)

        def first_valid_node(set_name):
            """Return the first independent node in the named set."""
            if set_name not in RVE_sets:
                return
            for node in RVE_sets[set_name].nodes:
                if (node.label not in secondary_labels
                        and node.label not in embedded_labels):
                    return node
            return

        # User-defined anchor
        if fixed_node is not None:
            set_name = (fixed_node if str(fixed_node).startswith("NSET_")
                        else "NSET_" + str(fixed_node))
            if set_name not in RVE_sets:
                raise RuntimeError(
                    "The node set '{}' selected for rigid-body suppression "
                    "does not exist.".format(set_name))
            node = first_valid_node(set_name)
            if node is None:
                raise RuntimeError(
                    "The node set '{}' selected for rigid-body suppression "
                    "does not contain an independent node that is free to be "
                    "constrained. Choose a non-embedded interior node or a "
                    "primary boundary node."
                    .format(set_name))
            return Region(nodes=MeshNodeArray(nodes=[node])), node

        # Automatic anchor selection
        pure_primary_groups = ["POINT_1", "EDGE_1", "EDGE_4", "EDGE_9",
                               "FACE_1", "FACE_3", "FACE_5"]
        for group in pure_primary_groups:
            node = first_valid_node("NSET_" + group)
            if node is not None:
                return Region(nodes=MeshNodeArray(nodes=[node])), node

        # Fallback: any boundary or inner node that is not secondary
        for set_name in ("NSET_BOUNDARY_NODES", "NSET_INNER_NODES",
                         "NSET_ALL_NODES"):
            node = first_valid_node(set_name)
            if node is not None:
                return Region(nodes=MeshNodeArray(nodes=[node])), node

        raise RuntimeError(
            "Could not find a suitable node to suppress rigid-body motion. "
            "Please provide one explicitly through the optional 'fixed_node' "
            "argument.")

    def assign_solid_section(self, element_set_name, material_name):
        """
        Assign a solid material section to an element set.

        Parameters
        ----------
        element_set_name : str
            The name of the element set that is assigned the section.
        material_name : str
            The name of the material that is associated with the
            section. The material should exist already.

        Returns
        -------
        None

        """

        self.model.HomogeneousSolidSection(name=element_set_name,
                                           material=material_name)
        self._assign_general_section(element_set_name, element_set_name)

    def assign_shell_section(self, element_set_name, material_name, thickness,
                             num_integration_points=5, pre_integrate=False):
        """
        Assign a homogeneous shell section to an element set.

        By default, section integration is performed during the analysis so
        nonlinear constitutive behavior and history variables are retained.

        Parameters
        ----------
        element_set_name : str
            The name of the element set that is assigned the section.
        material_name : str
            The name of the material that is associated with the
            section. The material should exist already.
        thickness : float
            The thickness of the shell.
        num_integration_points : int, optional
            The number of integration points through the thickness. The
            default is 5.
        pre_integrate : bool, optional
            If True, pre-integrate the section stiffness before the analysis;
            otherwise integrate through the thickness during the analysis.

        Returns
        -------
        None

        """

        integration = ON if pre_integrate else OFF
        self.model.HomogeneousShellSection(name=element_set_name,
                                           material=material_name,
                                           thickness=thickness,
                                           preIntegrate=integration,
                                           numIntPts=num_integration_points)
        self._assign_general_section(element_set_name, element_set_name)

    def assign_beam_section(self, element_set_name, material_name, radius,
                            pre_integrate=False):
        """
        Assign a beam section with circular profile to an element set.

        Beam orientations are assigned automatically. By default, section
        integration is performed during the analysis so nonlinear
        constitutive behavior and history variables are retained.

        Parameters
        ----------
        element_set_name : str
            The name of the element set that is assigned the section.
        material_name : str
            The name of the material that is associated with the
            section. The material should exist already.
        radius : float
            The radius of the circular cross-sectional area.
        pre_integrate : bool, optional
            If True, pre-integrate the section stiffness before the analysis;
            otherwise integrate over the cross-section during the analysis.

        Returns
        -------
        None

        """

        integration = BEFORE_ANALYSIS if pre_integrate else DURING_ANALYSIS
        self.model.CircularProfile(name=element_set_name, r=radius)
        self.model.BeamSection(name=element_set_name,
                               integration=integration,
                               profile=element_set_name,
                               material=material_name)
        self._assign_general_section(element_set_name, element_set_name)
        self._assign_beam_section_orientation(element_set_name)
        self._section_radius[element_set_name] = radius

    def assign_truss_section(self, element_set_name, material_name, area):
        """
        Assign a truss section to an element set.

        Parameters
        ----------
        element_set_name : str
            The name of the element set that is assigned the section.
        material_name : str
            The name of the material that is associated with the
            section. The material should exist already.
        area : float
            The cross-sectional area of the truss.

        Returns
        -------
        None

        """

        self.model.TrussSection(name=element_set_name, material=material_name,
                                area=area)
        self._assign_general_section(element_set_name, element_set_name)
        self._section_radius[element_set_name] = (area/np.pi)**0.5

    def _assign_general_section(self, section_name, element_set_name):
        """Assign an existing material section to an element set."""

        if element_set_name not in self.part.sets:
            element_set_name = "ELSET_{}".format(element_set_name)
        if element_set_name not in self.part.sets:
            raise RuntimeError("The element set was not found.")

        region = self.part.sets[element_set_name]
        self.part.SectionAssignment(region=region, sectionName=section_name)

    def _assign_beam_section_orientation(self, element_set_name):
        """Assign a cross-section direction orthogonal to each beam axis."""

        if element_set_name not in self.part.sets:
            element_set_name = "ELSET_{}".format(element_set_name)
        if element_set_name not in self.part.sets:
            raise RuntimeError("The element set was not found.")

        # For 2D beam models the first cross-section axis n1 must be the
        # out-of-plane direction
        if self.model_dimension == TWO_D_PLANAR:
            region = self.part.sets[element_set_name]
            self.part.assignBeamSectionOrientation(region=region,
                                                   method=N1_COSINES,
                                                   n1=(0.0, 0.0, -1.0))
            return

        # For 3D beam models the first cross-section axis n1 must be
        # perpendicular to the beam axis
        for element in self.part.sets[element_set_name].elements:
            # Compute unit tangent vector, i.e. the beam axis direction
            conn = element.connectivity
            node_start = self.part.nodes[conn[0]]
            node_end = self.part.nodes[conn[-1]]
            coords_start = np.array(node_start.coordinates)
            coords_end = np.array(node_end.coordinates)
            tangent = coords_end - coords_start
            tangent /= np.linalg.norm(tangent)

            # Pick a reference axis that is not parallel to the beam
            # axis and orthogonalize it against the tangent
            reference = np.array([1.0, 0.0, 0.0])
            if abs(np.dot(reference, tangent)) > 0.9:
                reference = np.array([0.0, 1.0, 0.0])
            normal = reference - np.dot(reference, tangent)*tangent
            normal /= np.linalg.norm(normal)

            assert abs(np.dot(tangent, normal)) < 1e-8

            # Assign orientation to the element
            region = Region(elements=MeshElementArray(elements=[element]))
            self.part.assignBeamSectionOrientation(region=region,
                                                   method=N1_COSINES,
                                                   n1=normal)

    def run_analysis(self, field_outputs=None, history_outputs=None,
                     **kwargs):
        """
        Create the job and run the finite element analysis in Abaqus.

        Parameters
        ----------
        field_outputs : iterable of str, optional
            The field outputs that shall be computed. The element nodal
            forces 'NFORC', nodal translations 'U', and nodal rotations
            'UR' needed by the homogenization post-processing are always
            requested in addition. The default is None.
        history_outputs : iterable of str, optional
            The history outputs that shall be computed. The default is None.
        **kwargs : object, optional
            Additional arguments that will be passed to the Abaqus job.
            See the documentation of the Abaqus Job class for details.

        Returns
        -------
        None

        """

        # Remove defaults
        if "F-Output-1" in self.model.fieldOutputRequests:
            del self.model.fieldOutputRequests["F-Output-1"]
        if "H-Output-1" in self.model.historyOutputRequests:
            del self.model.historyOutputRequests["H-Output-1"]

        # Do not modify containers supplied by the caller.
        field_outputs = [] if field_outputs is None else list(field_outputs)
        history_outputs = ([] if history_outputs is None
                           else list(history_outputs))

        # Request at least the outputs needed by average_stress and
        # average_power_increment.
        for output in ("NFORC", "U", "UR"):
            if output not in field_outputs:
                field_outputs.append(output)

        # Request the field outputs
        if len(field_outputs) > 0:
            for step in self.model.steps.keys():
                if step == "Initial":
                    continue
                name = "FIELD_OUTPUT_" + step
                if name in self.model.fieldOutputRequests:
                    del self.model.fieldOutputRequests[name]
                self.model.FieldOutputRequest(
                    name=name, createStepName=step, variables=field_outputs)

        # Request the history outputs
        if len(history_outputs) > 0:
            for step in self.model.steps.keys():
                if step == "Initial":
                    continue
                name = "HISTORY_OUTPUT_" + step
                if name in self.model.historyOutputRequests:
                    del self.model.historyOutputRequests[name]
                self.model.HistoryOutputRequest(
                    name=name, createStepName=step, variables=history_outputs)

        # Create and run the job. Repeated analyses use distinct names because
        # immediate reuse can leave solver resources associated with the old
        # job name, particularly in Abaqus Learning Edition.
        self.close_odb_file()
        self._job_count = getattr(self, "_job_count", 0) + 1
        job_name = "JOB_" + self.model_name
        if self._job_count > 1:
            job_name += "_{}".format(self._job_count)
        self._job_name = job_name
        if job_name in mdb.jobs.keys():
            del mdb.jobs[job_name]
        mdb.Job(name=job_name, model=self.model_name, **kwargs)
        mdb.jobs[job_name].submit()
        mdb.jobs[job_name].waitForCompletion()
        status = mdb.jobs[job_name].status
        if status is not None and status != COMPLETED:
            raise RuntimeError(
                "Abaqus job '{}' ended with status {}. See the job files "
                "for details.".format(job_name, status))

    def open_odb_file(self):
        """
        Open the output database file containing the simulation results.

        Post-processing methods reuse an output database opened through this
        method and leave it open when they return. If they open the database
        themselves, they also close it before returning.

        Returns
        -------
        None

        """

        odb_file = self.job_name + ".odb"
        if not os.path.isfile(odb_file):
            raise IOError(
                "The file {} does not exist".format(odb_file))

        if odb_file not in session.odbs.keys():
            session.openOdb(name=odb_file)
        self._odb = session.odbs[odb_file]

    def close_odb_file(self):
        """
        Close the output database file.

        Returns
        -------
        None

        """

        if self._odb is not None:
            self._odb.close()
        self._odb = None

    @staticmethod
    def _vtk_set_name(name, prefix):
        """Return an Abaqus set name with the required prefix."""

        return name if name.startswith(prefix) else prefix + name

    def _create_vtk_reader(self, field_outputs, node_sets, element_sets):
        """Create and configure the Paraqus output-database reader."""

        odb_file = self.job_name + ".odb"
        reader = OdbReader(odb_path=odb_file,
                           model_name=self.model_name,
                           instance_names=[self.model_name])
        for output in field_outputs:
            if type(output) is not str and len(output) == 2:
                reader.add_field_export_request(
                    output[0], invariant=output[1])
            else:
                reader.add_field_export_request(output)

        for name in element_sets:
            reader.add_set_export_request(
                self._vtk_set_name(name, "ELSET_"), set_type="ELEMENTS",
                instance_name=self.model_name)
        for name in node_sets:
            reader.add_set_export_request(
                self._vtk_set_name(name, "NSET_"), set_type="NODES",
                instance_name=self.model_name)
        return reader

    def _add_vtk_radius_field(self, instance_model, strut_radius):
        """Add the stored beam and truss radii to a Paraqus model."""

        if len(strut_radius) == 0:
            return

        element_radius = {}
        for name, value in strut_radius.items():
            set_name = self._vtk_set_name(name, "ELSET_")
            group = self.model_name + "." + set_name
            for tag in instance_model.elements.groups[group]:
                element_radius[tag] = value
        instance_model.add_field("RADIUS", list(element_radius.keys()),
                                 list(element_radius.values()),
                                 "ELEMENTS", "SCALAR")

    def _select_vtk_frames(self, reader, step_name=None, frame=None):
        """Return the requested ``(step_name, frame_index)`` pairs."""

        if frame is not None and (isinstance(frame, bool)
                                  or not isinstance(frame, Integral)):
            raise TypeError("The VTK frame index must be an integer.")

        available_steps = list(self.odb.steps.keys())
        if step_name is None:
            step_names = [name for name in available_steps
                          if name != "Initial"]
        elif step_name not in available_steps:
            raise ValueError(
                "The output database does not contain step '{}'."
                .format(step_name))
        else:
            step_names = [step_name]

        selection = []
        for name in step_names:
            available_frames = reader.get_frame_indices(name)
            if frame is None:
                selection.extend((name, index)
                                 for index in available_frames)
                continue

            frame_index = frame
            if frame_index < 0:
                frame_index += reader.get_number_of_frames(name)
            if frame_index not in available_frames:
                raise ValueError(
                    "Frame {} of step '{}' is unavailable or does not "
                    "contain all requested fields."
                    .format(frame, name))
            selection.append((name, frame_index))

        return selection

    def _write_vtk_frames(self, path, clear, reader, strut_radius,
                          step_name=None, frame=None):
        """Write the requested output frames through Paraqus."""

        vtu_writer = BinaryWriter(path, clear_output_dir=clear)
        with CollectionWriter(vtu_writer, self.model_name) as writer:
            selection = self._select_vtk_frames(
                reader, step_name=step_name, frame=frame)
            for current_step, frame_index in selection:
                instances = reader.read_instances(
                    step_name=current_step, frame_index=frame_index)
                instance_model = list(instances)[0]
                self._add_vtk_radius_field(instance_model, strut_radius)
                writer.write(instance_model)

    def export_to_vtk(self, path, field_outputs=None, node_sets=None,
                      element_sets=None, radius=None, clear=False,
                      step_name=None, frame=None):
        """
        Export the simulation results to VTK file format.

        Paraqus is used to read the Abaqus output database and write the VTK
        files.

        Parameters
        ----------
        path : str or os.PathLike
            Path to the folder in which the files will be stored.
        field_outputs : iterable of str or (str, str), optional
            Fields to export. An invariant is requested with a
            ``(field, invariant)`` tuple, e.g. ``("S", "MISES")``.
            The default is None.
        node_sets : iterable of str, optional
            The node sets that shall be exported. The default is None.
        element_sets : iterable of str, optional
            The element sets that shall be exported. The default is None.
        radius : mapping of str to float, optional
            A map from element set names to the strut radius. It overrides
            or supplements the radii recorded during section assignment.
            The default is None, i.e. only the recorded radii are used.
        clear : bool, optional
            If True, the output folder will be cleared if it already
            exists. The default is False.
        step_name : str, optional
            Name of the analysis step to export. If None, all non-initial
            steps are exported. The default is None.
        frame : int, optional
            Zero-based frame index to export. Negative indices count from
            the end of each step. If None, every frame containing all
            requested fields is exported. When ``step_name`` is None, the
            selected index is applied to every non-initial step. The default
            is None.

        Returns
        -------
        None

        """

        # Normalize the optional arguments (avoid mutable defaults)
        field_outputs = list(field_outputs) if field_outputs else []
        node_sets = list(node_sets) if node_sets else []
        element_sets = list(element_sets) if element_sets else []

        # Strut radii for the radius cell field
        strut_radius = dict(self._section_radius)
        if radius is not None:
            strut_radius.update(radius)

        # The strut element sets must be exported for the radius field
        for set_name in strut_radius.keys():
            if set_name not in element_sets:
                element_sets.append(set_name)

        if not os.path.isdir(path):
            os.makedirs(path)

        close_odb = self.odb is None
        self.open_odb_file()
        try:
            reader = self._create_vtk_reader(
                field_outputs, node_sets, element_sets)
            self._write_vtk_frames(
                path, clear, reader, strut_radius,
                step_name=step_name, frame=frame)
        finally:
            if close_odb:
                self.close_odb_file()

    def get_number_of_frames(self, step_name):
        """
        Return the number of frames for a given step.

        Parameters
        ----------
        step_name : str
            The name of the step.
        Returns
        -------
        n_frames : int
            The number of frames available in the given step.

        """

        close_odb = self.odb is None
        self.open_odb_file()
        try:
            return len(self.odb.steps[step_name].frames)
        finally:
            if close_odb:
                self.close_odb_file()

    def get_frame_time(self, step_name, frame=-1):
        """
        Return the total analysis time of an output frame.

        Abaqus stores a frame value relative to the beginning of its step.
        The total time is therefore the sum of the time accumulated before
        the step and the frame value within the step.

        Parameters
        ----------
        step_name : str
            Name of the analysis step.
        frame : int, optional
            Frame index. The default is -1.

        Returns
        -------
        float
            Total analysis time at the requested frame.

        """

        close_odb = self.odb is None
        self.open_odb_file()
        try:
            step = self.odb.steps[step_name]
            if step.domain != TIME:
                raise ValueError(
                    "Frame time is only defined for time-domain steps.")
            return float(step.totalTime + step.frames[frame].frameValue)
        finally:
            if close_odb:
                self.close_odb_file()

    def average_stress(self, step_name, frame=-1):
        r"""
        Compute the effective first Piola--Kirchhoff stress.

        With the force-on-node convention used by Abaqus ``NFORC``, the
        element nodal-force virial of the retained reduced model is

        .. math::
            \bar{\boldsymbol{P}}=-\frac{1}{V}
            \sum_{e\in\mathcal{S}_{h}^{\mathrm{red}}}
            \sum_{j\in e}\boldsymbol{f}^{(e)}_{j}
            \otimes\boldsymbol{X}_{j}.

        In the geometrically linear regime,
        :math:`\bar{\boldsymbol{P}}\approx\bar{\boldsymbol{\sigma}}`, so the
        same expression returns the effective Cauchy stress to first order.

        Parameters
        ----------
        step_name : str
            Name of the step.
        frame : int, optional
            The frame index. The default is -1.

        Returns
        -------
        stress : ndarray
            Effective first Piola--Kirchhoff stress
            :math:`\bar{\boldsymbol{P}}` in full tensor form, with shape
            ``(2, 2)`` or ``(3, 3)``.

        """

        close_odb = self.odb is None
        self.open_odb_file()
        try:
            dim = 3 if self.model_dimension == THREE_D else 2
            volume = self.get_volume()
            output_frame = self.odb.steps[step_name].frames[frame]

            # NFORC is stored as scalar components
            field_outputs = output_frame.fieldOutputs
            if "NFORC1" not in field_outputs.keys():
                raise RuntimeError("Field output 'NFORC' is not available.")

            odb_instance = self.odb.rootAssembly.instances[self.model_name]
            node_coords = {n.label: np.array(n.coordinates[0:dim])
                           for n in odb_instance.nodes}

            # Assemble the internal nodal force vectors per (element, node)
            forces = {}
            for i in range(dim):
                for f in field_outputs["NFORC{}".format(i + 1)].values:
                    if f.elementLabel is None:
                        continue
                    key = (f.elementLabel, f.nodeLabel)
                    if key not in forces:
                        forces[key] = np.zeros(dim)
                    forces[key][i] = f.data

            # Virial sum over the retained RVE elements
            stress = np.zeros((dim, dim), dtype=float)
            for (_, node_label), force in forces.items():
                stress += np.outer(force, node_coords[node_label])
            stress *= -1.0/volume
            return stress
        finally:
            if close_odb:
                self.close_odb_file()

    @staticmethod
    def _hat_operator(vector):
        """Return the skew tensor associated with an axial vector."""

        return np.array([[0.0, -vector[2], vector[1]],
                         [vector[2], 0.0, -vector[0]],
                         [-vector[1], vector[0], 0.0]])

    @staticmethod
    def _left_jacobian(rotation):
        """Evaluate the SO(3) left Jacobian, using its series near zero."""

        angle = np.linalg.norm(rotation)
        rotation_hat = AbaqusRVE._hat_operator(rotation)
        if angle < 1.0e-8:
            return (np.eye(3) + 0.5*rotation_hat
                    + rotation_hat.dot(rotation_hat)/6.0)
        return (np.eye(3)
                + (1.0 - np.cos(angle))/angle**2*rotation_hat
                + (angle - np.sin(angle))/angle**3
                * rotation_hat.dot(rotation_hat))

    @staticmethod
    def _odb_vector_field(outputs, name, instance):
        """Return an ODB nodal vector field keyed by node label."""

        if name not in outputs:
            return {}
        field = outputs[name].getSubset(region=instance)
        return {value.nodeLabel: np.atleast_1d(
                    np.asarray(value.data, dtype=float))
                for value in field.values}

    def _odb_rotation_field(self, outputs, instance, dim):
        """Return ODB nodal rotations as three-component vectors."""

        if "UR" in outputs:
            rotations = self._odb_vector_field(outputs, "UR", instance)
            result = {}
            for node, value in rotations.items():
                if len(value) < 3:
                    raise RuntimeError(
                        "Rotational output 'UR' has fewer than three "
                        "components.")
                result[node] = np.asarray(value[:3], dtype=float)
            return result

        component_names = (("UR1", "UR2", "UR3") if dim == 3
                           else ("UR3",))
        result = {}
        for local_index, name in enumerate(component_names):
            if name not in outputs:
                continue
            field = outputs[name].getSubset(region=instance)
            component = local_index if dim == 3 else 2
            for value in field.values:
                result.setdefault(value.nodeLabel, np.zeros(3))
                result[value.nodeLabel][component] = float(value.data)
        return result

    @staticmethod
    def _consecutive_frames(step, frame):
        """Return the selected output frame and its predecessor."""

        n_frames = len(step.frames)
        frame_index = frame if frame >= 0 else n_frames + frame
        if frame_index <= 0 or frame_index >= n_frames:
            raise IndexError(
                "The selected frame must have a preceding output frame.")
        return step.frames[frame_index], step.frames[frame_index - 1]

    @staticmethod
    def _translation_increments(current, previous, dim):
        """Return nodal translation increments between two frames."""

        increments = {}
        for node, value in current.items():
            if node not in previous:
                raise RuntimeError(
                    "Nodal translation output is inconsistent between "
                    "the consecutive frames.")
            increments[node] = value[:dim] - previous[node][:dim]
        return increments

    def _rotation_increments(self, current, previous):
        """Map rotation-vector increments to current spatial increments."""

        increments = {}
        for node, rotation in current.items():
            if node not in previous:
                raise RuntimeError(
                    "Nodal rotation output is inconsistent between the "
                    "consecutive frames.")
            increment = rotation - previous[node]
            increments[node] = self._left_jacobian(rotation).dot(increment)
        return increments

    @staticmethod
    def _has_nodal_moments(outputs, instance, dim):
        """Return whether the current frame contains nodal moments."""

        names = (["NFORC4", "NFORC5", "NFORC6"] if dim == 3
                 else ["NFORC6"])
        return any(name in outputs
                   and len(outputs[name].getSubset(
                       region=instance).values) > 0
                   for name in names)

    @staticmethod
    def _translation_incremental_work(outputs, instance, increments, dim):
        """Assemble the translational contribution to incremental work."""

        work = 0.0
        for component in range(dim):
            name = "NFORC{}".format(component + 1)
            if name not in outputs:
                raise RuntimeError(
                    "Field output '{}' is required.".format(name))
            field = outputs[name].getSubset(region=instance)
            for value in field.values:
                if value.elementLabel is None:
                    continue
                increment = increments.get(value.nodeLabel)
                if increment is not None:
                    work -= float(value.data)*increment[component]
        return work

    @staticmethod
    def _rotation_incremental_work(outputs, instance, increments, dim):
        """Assemble the moment--rotation contribution to incremental work."""

        work = 0.0
        components = ((3, 4, 5) if dim == 3 else (5,))
        for component in components:
            name = "NFORC{}".format(component + 1)
            if name not in outputs:
                continue
            field = outputs[name].getSubset(region=instance)
            rotation_component = component - 3 if dim == 3 else 2
            for value in field.values:
                if value.elementLabel is None:
                    continue
                rotation = increments.get(value.nodeLabel)
                if rotation is not None:
                    work -= float(value.data)*rotation[rotation_component]
        return work

    def average_power_increment(self, step_name, frame=-1):
        r"""
        Compute the microscopic incremental-work density.

        The finite increments between the selected frame and the immediately
        preceding frame are inserted as a test direction in the discrete
        Hill--Mandel relation. With current-state element forces and moments,

        .. math::
            \Delta w_\mathrm{micro} = -\frac{1}{V}
            \sum_{e\in\mathcal{S}_{h}^{\mathrm{red}}}\sum_{j\in e}
            \left(\boldsymbol{f}^{(e)}_j\cdot\Delta\boldsymbol{u}_j
            +\boldsymbol{m}^{(e)}_j\cdot
            \boldsymbol{\eta}^{\Delta}_j\right),

        where :math:`\boldsymbol{\eta}^{\Delta}_j =
        \boldsymbol{J}_\mathrm{l}(\boldsymbol{\phi}_j)
        \Delta\boldsymbol{\phi}_j`. In two-dimensional models, only the
        out-of-plane rotational component is considered. Nodes without
        rotational degrees of freedom contribute only through their
        translations.

        The analysis must contain the field outputs ``NFORC``, ``U`` and,
        when rotational degrees of freedom are present, ``UR``. These are
        requested by default by :meth:`run_analysis`.

        Parameters
        ----------
        step_name : str
            Name of the analysis step.
        frame : int, optional
            Current frame index. The immediately preceding frame supplies
            the finite increment. The default is -1.

        Returns
        -------
        float
            Microscopic incremental-work density
            :math:`\Delta w_\mathrm{micro}`.

        """

        close_odb = self.odb is None
        self.open_odb_file()
        try:
            step = self.odb.steps[step_name]
            current, previous = self._consecutive_frames(step, frame)
            current_outputs = current.fieldOutputs
            previous_outputs = previous.fieldOutputs
            if "U" not in current_outputs or "U" not in previous_outputs:
                raise RuntimeError(
                    "Field output 'U' is required in both consecutive "
                    "frames.")
            if "NFORC1" not in current_outputs:
                raise RuntimeError(
                    "Field output 'NFORC' is required in the current frame.")

            dim = 3 if self.model_dimension == THREE_D else 2
            odb_instance = self.odb.rootAssembly.instances[self.model_name]

            translations_current = self._odb_vector_field(
                current_outputs, "U", odb_instance)
            translations_previous = self._odb_vector_field(
                previous_outputs, "U", odb_instance)
            rotations_current = self._odb_rotation_field(
                current_outputs, odb_instance, dim)
            rotations_previous = self._odb_rotation_field(
                previous_outputs, odb_instance, dim)

            if (self._has_nodal_moments(current_outputs, odb_instance, dim)
                    and (not rotations_current or not rotations_previous)):
                raise RuntimeError(
                    "Field output 'UR' is required in both consecutive "
                    "frames when nodal moments are present.")

            translation_increment = self._translation_increments(
                translations_current, translations_previous, dim)
            rotation_increment = self._rotation_increments(
                rotations_current, rotations_previous)

            incremental_work = self._translation_incremental_work(
                current_outputs, odb_instance, translation_increment, dim)
            incremental_work += self._rotation_incremental_work(
                current_outputs, odb_instance, rotation_increment, dim)

            return incremental_work/self.get_volume()
        finally:
            if close_odb:
                self.close_odb_file()


class LoadHistory(object):
    r"""
    Tensor-valued macroscopic loading history.

    For the boundary-condition methods, the returned tensor is the
    macroscopic displacement gradient
    :math:`\bar{\boldsymbol{H}}(\tau)=
    \bar{\boldsymbol{F}}(\tau)-\boldsymbol{I}` sampled at total analysis
    times representing the pseudo-time :math:`\tau`.

    Parameters
    ----------
    time : one-dimensional ArrayLike
        Finite, strictly increasing total analysis-time values at which the
        history is prescribed. At least one value is required.
    history : callable
        Function mapping time to a tensor represented by a ``(3, 3)`` array.

    Attributes
    ----------
    time_data : ndarray
        Stored total analysis-time values.
    load_data : ndarray
        Tensor values evaluated at the stored times.

    """

    def __init__(
            self, time, history):
        self._time = np.asarray(time, dtype=float)
        if self._time.ndim != 1:
            raise ValueError("The time data must be one-dimensional.")
        if self._time.size == 0:
            raise ValueError("The time data must contain at least one value.")
        if not np.all(np.isfinite(self._time)):
            raise ValueError("The time data must contain only finite values.")
        if np.any(np.diff(self._time) <= 0.0):
            raise ValueError("The time data must be strictly increasing.")
        if not callable(history):
            raise TypeError("The history must be callable.")
        self._history = history
        self._index = 0

        # Test return value of history
        self[self._index]

    @property
    def time_data(self):
        """
        Total analysis-time values of the load history.

        Returns
        -------
        ndarray
            The time data points.

        """
        return self._time

    @property
    def load_data(self):
        """
        Tensor values evaluated at all stored time values.

        Returns
        -------
        data : ndarray
            The load data points.

        """
        data = []
        for _, v in self:
            data.append(v)
        data = np.array(data)
        return data

    def __iter__(self):
        self._index = 0
        return self

    def __next__(self):
        if self._index >= len(self.time_data):
            raise StopIteration

        time, load = self[self._index]
        self._index += 1

        return time, load

    def next(self):
        """
        Return the next time--load pair for legacy Abaqus with Python 2.7.

        Returns
        -------
        time : float
            Current total analysis time.
        load : ndarray
            Macroscopic displacement gradient at the current time.

        """

        time, load = self.__next__()
        return time, load

    def __getitem__(self, i):
        time = float(self.time_data[i])
        load = np.array(self(time))

        if load.shape != (3, 3):
            raise ValueError(
                "The history function must return the load in full 3D tensor "
                "notation as a (3, 3) array.")

        return time, load

    def __call__(self, time):
        """
        Evaluate the prescribed loading history.

        Parameters
        ----------
        time : float
            Total analysis time at which the history is evaluated.

        Returns
        -------
        load : ArrayLike
            Macroscopic displacement gradient returned by the history
            function.

        """
        return self._history(time)

    def __str__(self):
        descr = "LoadHistory with {} data points\n".format(len(self.time_data))
        return descr

    def __repr__(self):
        return self.__str__()

    def as_table(self):
        """
        Get all pairs of time and load at once as a table.

        Returns
        -------
        data : tuple of 2-tuples
            The ``(time, load)`` pairs, in which each time is a float and
            each load is an ndarray.

        """

        data = []
        for t, v in self:
            data.append((t, np.array(v)))
        data = tuple(data)
        return data
