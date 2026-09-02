"""
Module containing the AbaqusRVE class that provides useful methods for
computational homogenization with Abaqus.

"""

import os
import numpy as np

from abqhom.constants import (PERIODIC_PAIRS, ELEMENTS_ROTATIONAL_DOFS,
                              PARAQUS_INSTALLED, ELEMENTS_WARPING_DOFS)

from abaqus import Mdb, mdb, session
from regionToolset import Region
from mesh import MeshNodeArray, MeshElementArray
from abaqusConstants import (THREE_D, TWO_D_PLANAR, DEFORMABLE_BODY, ON, OFF,
                             UNSET, CARTESIAN, N1_COSINES, DURING_ANALYSIS,
                             BEFORE_ANALYSIS, TOTAL)

if PARAQUS_INSTALLED:
    from paraqus.abaqus import OdbReader
    from paraqus.writers import BinaryWriter, CollectionWriter


class AbaqusRVE(object):
    """
    An RVE representation to set up a homogenization problem in Abaqus.

    Parameters
    ----------
    gmsh_RVE : GmshRVE
        The GmshRVE containing the mesh information of the RVE.
    custom_types : dict[int, str], optional
        Map from Gmsh element types to Abaqus element types. If this
        is an empty dict a default map will be used. The default is
        {}.

    Attributes
    ----------
    model_name : str
        The name of the model.
    periodic_mesh : bool
        A flag to mark whether the model has a mesh with periodicity
        constraints.
    model_dimension : SymbolicConstant
        The model dimension. The value can be THREE_D or TWO_D_PLANAR.
    model : Model
        The internal Abaqus representation of the model.
    part : Part
        The internal Abaqus representation of the RVE part.
    instance : PartInstance
        The internal representation of the RVE part in the assembly.
    odb : Odb
        The Abaqus output database containing the analysis results.

    Example
    -------
    >>> from abqhom.gmshmodel import GmshRVE
    >>> from abqhom.abqmodel import AbaqusRVE, LoadHistory
    >>>
    >>> # Assume a model with the name GMSH_MODEL exists in Gmsh
    >>> gmsh_RVE = GmshRVE("GMSH_MODEL")
    >>> gmsh_RVE.add_periodic_mesh_constraints()
    >>> gmsh_RVE.generate_mesh()
    >>>
    >>> # Strain in x-direction
    >>> def strain_path(time):
    >>>     strain = np.zeros((3, 3))
    >>>     strain[0, 0] = 0.1*time
    >>>     return strain
    >>>
    >>> # Create the load history
    >>> time = np.linspace(0.0, 1.0, 50)
    >>> history = LoadHistory(time, strain_path)
    >>>
    >>> # Create an AbaqusRVE from the GmshRVE (the input file is written
    >>> # and imported automatically)
    >>> unit_cell = AbaqusRVE(gmsh_RVE)
    >>>
    >>> # Assign a material and a solid section to an element set
    >>> material = unit_cell.model.Material(name="Steel")
    >>> material.Elastic(table=[[210.0, 0.3]])
    >>> unit_cell.assign_solid_section("ALL_ELEMENTS", "Steel")
    >>>
    >>> # Add a step, apply periodic boundary conditions and run
    >>> unit_cell.model.StaticStep(name="Load", previous="Initial")
    >>> unit_cell.apply_periodic_BC("Load", history)
    >>> unit_cell.run_analysis(field_outputs=["U", "S"])
    >>>
    >>> # Compute the effective stress and export the results to VTK
    >>> stress = unit_cell.average_stress("Load")
    >>> unit_cell.export_to_vtk(vtk_path, field_outputs=["U", ("S", "MISES")])

    """

    def __init__(self, gmsh_RVE, custom_types={}):
        if gmsh_RVE.model_dimension not in {2, 3}:
            raise RuntimeError("Only 2d and 3d models are supported.")
        dim = THREE_D if gmsh_RVE.model_dimension == 3 else TWO_D_PLANAR

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
        A flag to mark whether the mesh has periodicity constraints.

        Returns
        -------
        bool
            True if the mesh is periodic, False otherwise.

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
        Odb
            The output database.

        """
        return self._odb

    def _add_node_sets(self, node_sets):
        """
        Add node sets to the Abaqus model.

        Parameters
        ----------
        node_sets : dict[str, Iterable]
            The node sets with the set names as keys and the set nodes
            as values.

        Returns
        -------
        None.

        """

        for set_name, set_items in node_sets.items():
            set_name = "NSET_" + str(set_name)
            set_items = [int(i) for i in set_items]

            # Sets must be unsorted to not break ordering for periodic BCs
            if len(set_items) > 0:
                self.part.SetFromNodeLabels(name=set_name,
                                            nodeLabels=set_items,
                                            unsorted=True)

    def _add_element_sets(self, element_sets):
        """
        Add element sets to the Abaqus model.

        Parameters
        ----------
        element_sets : dict[str, Iterable]
            The element sets with the set names as keys and the set
            elements as values.

        Returns
        -------
        None.

        """

        for set_name, set_items in element_sets.items():
            set_name = "ELSET_" + str(set_name)
            set_items = [int(i) for i in set_items]

            if len(set_items) > 0:
                self.part.SetFromElementLabels(name=set_name,
                                               elementLabels=set_items)

    def get_bounding_box(self):
        """
        Return the geometric bounding box of the RVE.

        The bounding box is the one stored at construction (computed from
        the Gmsh CAD points), not recomputed from the mesh.

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

    def _collect_element_types(self):
        """
        Create a set containing the names of all occuring element types.

        Returns
        -------
        element_types : set[str]
            The string representations of all element types in the mesh.

        """

        all_elements = self.part.sets["ELSET_ALL_ELEMENTS"].elements
        element_types = {str(element.type) for element in all_elements}
        return element_types

    def _add_dirichlet_from_displacement_history(self, time_history, u_history,
                                                 region, step_name, bc_name,
                                                 amplitude_name):
        """
        Add history-dependent Dirichlet boundary conditions.

        Parameters
        ----------
        time_history : ndarray
            The time data points of the displacement history.
        u_history : ndarray
            The displacement data points of the displacement history.
        region : Region
            The node region for applying the boundary conditions.
        step_name : str
            The name of the step in which the boundary conditions will
            be initialized.
        bc_name : str
            The name of the boundary condition.
        amplitude_name : str
            The name of the displacement history amplitude.

        Returns
        -------
        None.

        """

        for i in range(3):
            # Ignore z component for 2d models
            if i == 2 and self.model_dimension != THREE_D:
                break

            # Add amplitude for the load history
            amplitude_name_i = amplitude_name + "_{}".format(i + 1)
            data = np.vstack((time_history, u_history[:, i])).T
            data = tuple([(t, v) for t, v in data])
            self.model.TabularAmplitude(name=amplitude_name_i, timeSpan=TOTAL,
                                        data=data)

            # Add the boundary conditon
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

    def apply_affine_displacement_BC(self, initialization_step,
                                     strain_history):
        r"""
        Apply affine displacement boundary conditions.

        The boundary conditions will be transferred to all subsequent
        steps, i.e. the strain history must be formulated w.r.t. the
        total simulation time and not the step time.

        Affine displacement boundary conditions read

        .. math::
            \boldsymbol{u} =
            \boldsymbol{\mathcal{E}} \cdot \boldsymbol{X}

        whereby the prescribed macroscopic strain
        :math:`\boldsymbol{\mathcal{E}}`
        is the linearized strain tensor
        :math:`\boldsymbol{\mathcal{E}} = \boldsymbol{\varepsilon}`
        for geometrically linear analyses and the deformation gradient
        minus the identity tensor
        :math:`\boldsymbol{\mathcal{E}} = \boldsymbol{F} -
        \boldsymbol{I}`
        for geometrically nonlinear analyses.

        **Elements with rotational or warping degrees of freedom (beams,
        shells, open-section beams) are not supported by the affine
        displacement boundary conditions.**

        Parameters
        ----------
        initialization_step : str
            The name of the step in which the boundary conditions are
            initialized.
        strain_history : LoadHistory
            The load history containing time and strain data.

        Returns
        -------
        None.

        """

        if self._bc_initialized:
            raise RuntimeError(
                "Boundary conditions have already been initialized.")

        # Rotational/warping DOFs (beams, shells, open-section beams) are
        # not supported: the consistent boundary rotation depends on the
        # macroscopic rotation, which is only defined once a macroscopic
        # rotation/curvature model is fixed. Use periodic BCs instead.
        rotational_or_warping = ELEMENTS_ROTATIONAL_DOFS | ELEMENTS_WARPING_DOFS
        if len(self._collect_element_types() & rotational_or_warping) > 0:
            raise NotImplementedError(
                "Affine displacement boundary conditions are not supported "
                "for elements with rotational or warping degrees of freedom "
                "(beams, shells, open-section beams). Please use periodic "
                "boundary conditions instead.")

        time_history = strain_history.time_data
        strain_history = strain_history.load_data

        # Apply Dirichlet boundary conditions on all boundary nodes
        boundary_nodes = self.instance.sets["NSET_BOUNDARY_NODES"]
        for node in boundary_nodes.nodes:
            region = Region(nodes=MeshNodeArray(nodes=[node]))
            x = np.array(node.coordinates)
            u_history = np.dot(strain_history, x)

            amplitude_name = "DIRICHLET_{}".format(node.label)
            bc_name = "DIRICHLET_{}".format(node.label)
            self._add_dirichlet_from_displacement_history(
                time_history, u_history, region, initialization_step, bc_name,
                amplitude_name)

        self._bc_initialized = True

    def apply_periodic_BC(self, initialization_step, strain_history,
                          fixed_node=None):
        r"""
        Apply periodic displacement boundary conditions.

        The boundary conditions will be transferred to all subsequent
        steps, i.e. the strain history must be formulated w.r.t. the
        total simulation time and not the step time.

        Periodic displacement boundary conditions read

        .. math::
            \boldsymbol{u}^{+} - \boldsymbol{u}^{-} =
            \boldsymbol{\mathcal{E}} \cdot \left[
                \boldsymbol{X}^{+} - \boldsymbol{X}^{-} \right]

        whereby the prescribed macroscopic strain
        :math:`\boldsymbol{\mathcal{E}}`
        is the linearized strain tensor
        :math:`\boldsymbol{\mathcal{E}} = \boldsymbol{\varepsilon}`
        for geometrically linear analyses and the deformation gradient
        minus the identity tensor
        :math:`\boldsymbol{\mathcal{E}} = \boldsymbol{F} -
        \boldsymbol{I}`
        for geometrically nonlinear analyses.

        If the model contains rotational degrees of freedom, additional
        constraints of the form

        .. math::
            \boldsymbol{\theta}^{+} - \boldsymbol{\theta}^{-} =
            \boldsymbol{0}

        will be enforced, implying that the macroscopic curvature is
        assumed to be zero.

        Parameters
        ----------
        initialization_step : str
            The name of the step in which the boundary conditions are
            initialized.
        strain_history : LoadHistory
            The load history containing time and strain data.
        fixed_node : str, optional
            The name of a node set whose first admissible node is used
            as the anchor. If None, an anchor is selected automatically.
            The default is None.

        Returns
        -------
        None.

        """

        if self._bc_initialized:
            raise RuntimeError(
                "Boundary conditions have already been initialized.")

        time_history = strain_history.time_data
        strain_history = strain_history.load_data

        def create_reference_point(master, slave):
            """
            Create a reference point to connect master and slave nodes.

            Parameters
            ----------
            master : str
                The name of the master boundary group.
            slave : str
                The name of the slave boundary group.

            Returns
            -------
            None.

            """

            coords = [0.0, 0.0, 0.0]
            point_name = "REFERENCE_POINT_{}_{}".format(master, slave)
            part = self.model.Part(dimensionality=self.model_dimension,
                                   name=point_name, type=DEFORMABLE_BODY)
            point = self.model.parts[point_name].ReferencePoint(point=coords)
            self.model.parts[point_name].Set(
                name="REFERENCE_POINT",
                referencePoints=[self.model.parts[point_name]
                                 .referencePoints[point.id]])
            self.model.rootAssembly.Instance(name=point_name, part=part,
                                             dependent=ON)

            # Compute the +/- coordinates of the resepctive groups
            x_minus = RVE_sets["NSET_" + master].nodes[0].coordinates
            x_plus = RVE_sets["NSET_" + slave].nodes[0].coordinates

            # Apply the boundary conditions
            region = self.model.rootAssembly.instances[point_name].sets[
                "REFERENCE_POINT"]
            u_history = np.dot(strain_history,
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
        region, anchor_node = self._select_rigid_body_anchor(fixed_node)
        x = np.array(anchor_node.coordinates)
        u_history = np.dot(strain_history, x)
        amplitude_name = "DIRICHLET_ANCHOR_{}".format(anchor_node.label)
        bc_name = "DIRICHLET_ANCHOR_{}".format(anchor_node.label)
        self._add_dirichlet_from_displacement_history(
            time_history, u_history, region, initialization_step, bc_name,
            amplitude_name)

        # Add periodic coupling to edges and faces
        for slave, master in PERIODIC_PAIRS.items():
            if ("NSET_" + master not in RVE_sets
                    or "NSET_" + slave not in RVE_sets):
                continue

            # Create a reference point for the RHS of the equation
            create_reference_point(master, slave)

            for dof in range(1, dim + 1):
                base_name = self.model_name + ".NSET_"
                ref_set_name = ("REFERENCE_POINT_{}_{}.REFERENCE_POINT"
                                .format(master, slave))
                terms = [[1.0, base_name + slave, dof],
                         [-1.0, base_name + master, dof],
                         [-1.0, ref_set_name, dof]]

                self.model.Equation(name="PBC_{}_{}_{}".format(
                    master, slave, dof), terms=terms)

        # Apply constraints also to rotations and warping
        self._periodic_boundary_rotations()

        self._bc_initialized = True
        self._has_periodic_bc = True

    def _periodic_boundary_rotations(self):
        """
        Apply periodic boundary conditions to higher-order DOFs.

        Currently considered higher-order DOFs are rotational DOFs
        (beams) and warping DOFs (open section beams). Coupling the
        slave and master DOFs without a reference-point term implies a
        vanishing macroscopic curvature.

        Returns
        -------
        None.

        """

        # Collect all element types present in the model
        all_element_types = self._collect_element_types()

        # Collect the dofs that need to be constrained
        dofs = []
        if len(all_element_types & ELEMENTS_ROTATIONAL_DOFS) > 0:
            rotational_dofs = ([4, 5, 6] if self.model_dimension == THREE_D
                               else [6])
            dofs.extend(rotational_dofs)  # Rotations

        if len(all_element_types & ELEMENTS_WARPING_DOFS) > 0:
            dofs.append(7)  # Warping

        # Apply equation constraints to rotational degrees of freedom
        for slave, master in PERIODIC_PAIRS.items():
            if ("NSET_" + master not in self.instance.sets
                    or "NSET_" + slave not in self.instance.sets):
                continue

            for dof in dofs:
                base_name = self.model_name + ".NSET_"
                terms = [[1.0, base_name + slave, dof],
                         [-1.0, base_name + master, dof]]

                self.model.Equation(name="PBC_{}_{}_{}".format(
                    master, slave, dof), terms=terms)

    def _select_rigid_body_anchor(self, fixed_node=None):
        r"""
        Select a single node to suppress rigid-body motion.

        Under periodic boundary conditions the macroscopic strain fixes
        all relative boundary displacements, so only the rigid-body
        translations remain. They are removed by prescribing the affine
        displacement :math:`\boldsymbol{u} = \boldsymbol{\mathcal{E}}
        \cdot \boldsymbol{X}` at a single anchor node. The anchor must not
        be a dependent (slave) node of a periodic equation, otherwise
        Abaqus reports an over-constraint. Nodes of the "pure master"
        boundary groups (which never act as a slave) and interior nodes
        satisfy this requirement.

        Parameters
        ----------
        fixed_node : str, optional
            The name of a node set whose first admissible node is used
            as the anchor. If None, an anchor is selected automatically.
            The default is None.

        Returns
        -------
        region : Region
            A region containing exactly the single anchor node.
        node : MeshNode
            The anchor node.

        """

        RVE_sets = self.instance.sets

        # Collect the labels of all periodic slave nodes
        slave_labels = set()
        for slave in PERIODIC_PAIRS.keys():
            set_name = "NSET_" + slave
            if set_name in RVE_sets:
                slave_labels.update(n.label for n in RVE_sets[set_name].nodes)

        def first_valid_node(set_name):
            if set_name not in RVE_sets:
                return
            for node in RVE_sets[set_name].nodes:
                if node.label not in slave_labels:
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
                    "does not contain a node that is free to be constrained. "
                    "Choose an interior node or a master boundary node."
                    .format(set_name))
            return Region(nodes=MeshNodeArray(nodes=[node])), node

        # Automatic anchor selection
        pure_master_groups = ["POINT_1", "EDGE_1", "EDGE_4", "EDGE_9",
                              "FACE_1", "FACE_3", "FACE_5"]
        for group in pure_master_groups:
            node = first_valid_node("NSET_" + group)
            if node is not None:
                return Region(nodes=MeshNodeArray(nodes=[node])), node

        # Fallback: any boundary or inner node that is not a slave
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
        None.

        """

        self.model.HomogeneousSolidSection(name=element_set_name,
                                           material=material_name)
        self._assign_general_section(element_set_name, element_set_name)

    def assign_shell_section(self, element_set_name, material_name, thickness,
                             num_integration_points=5, pre_integrate=False):
        """
        Assign a homogeneous shell section to an element set.

        By default the section is integrated through the thickness during
        the analysis (``preIntegrate=OFF``) rather than pre-integrated, so
        that nonlinear material behaviour is captured correctly --
        analogous to the beam section integrated during the analysis.

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
            If False (the default) the section is integrated through the
            thickness during the analysis, supporting nonlinear material
            behaviour. If True the section stiffness is pre-integrated
            before the analysis, which is cheaper but restricted to
            linear material behaviour. The default is False.

        Returns
        -------
        None.

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

        The orientation of the beams will be created automatically. By
        default the section is integrated over the cross-section during
        the analysis, so that nonlinear material behaviour is captured.

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
            If False (the default) the section is integrated over the
            cross-section during the analysis, supporting nonlinear
            material behaviour. If True the section stiffness is
            pre-integrated before the analysis, which is cheaper but
            restricted to linear material behaviour. The default is False.

        Returns
        -------
        None.

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
            The cross-sectional area of the beam.

        Returns
        -------
        None.

        """

        self.model.TrussSection(name=element_set_name, material=material_name,
                                area=area)
        self._assign_general_section(element_set_name, element_set_name)
        self._section_radius[element_set_name] = (area/np.pi)**0.5

    def _assign_general_section(self, section_name, element_set_name):
        """
        Assign any material section to an element set.

        Parameters
        ----------
        section_name : str
            The name of the material section.
        element_set_name : str
            The name of the element set that is assigned the section.

        Returns
        -------
        None.

        """

        if element_set_name not in self.part.sets:
            element_set_name = "ELSET_{}".format(element_set_name)
        if element_set_name not in self.part.sets:
            raise RuntimeError("The element set was not found.")

        region = self.part.sets[element_set_name]
        self.part.SectionAssignment(region=region, sectionName=section_name)

    def _assign_beam_section_orientation(self, element_set_name):
        """
        Compute the correct tangent and normal vectors of beams.

        This method is only valid for beam elements with a circular
        profile.

        Parameters
        ----------
        element_set_name : str
            The name of the element set that is assigned a beam section.

        Returns
        -------
        None.

        """

        if element_set_name not in self.part.sets:
            element_set_name = "ELSET_{}".format(element_set_name)
        if element_set_name not in self.part.sets:
            raise RuntimeError("The element set was not found.")

        # For 2d beam models the first cross-section axis n1 must be the
        # out-of-plane direction
        if self.model_dimension == TWO_D_PLANAR:
            region = self.part.sets[element_set_name]
            self.part.assignBeamSectionOrientation(region=region,
                                                   method=N1_COSINES,
                                                   n1=(0.0, 0.0, -1.0))
            return

        # For 3d beam models the first cross-section axis n1 must be
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

    def run_analysis(self, field_outputs=[], history_outputs=[],
                     **kwargs):
        """
        Create the job and run the finite element analysis in Abaqus.

        Parameters
        ----------
        field_outputs : list[str], optional
            The field outputs that shall be computed. The element nodal
            forces 'NFORC' needed by average_stress are always requested
            in addition. The default is [].
        history_outputs : list[str], optional
            The history outputs that shall be computed. The default is
            [].
        **kwargs :
            Additional arguments that will be passed to the Abaqus job.
            See the documentation of the Abaqus Job class for details.

        Returns
        -------
        None.

        """

        # Remove defaults
        if "F-Output-1" in self.model.fieldOutputRequests:
            del self.model.fieldOutputRequests["F-Output-1"]
        if "H-Output-1" in self.model.historyOutputRequests:
            del self.model.historyOutputRequests["H-Output-1"]

        # Request at least the outputs needed by average_stress
        if "NFORC" not in field_outputs:
            field_outputs.append("NFORC")

        # Request the field outputs
        if len(field_outputs) > 0:
            for step in self.model.steps.keys():
                if step == "Initial":
                    continue
                self.model.FieldOutputRequest(
                    name="FIELD_OUTPUT_" + step, createStepName=step,
                    variables=field_outputs)
                    # timeInterval=EVERY_TIME_INCREMENT)

        # Request the history outputs
        if len(history_outputs) > 0:
            for step in self.model.steps.keys():
                if step == "Initial":
                    continue
                self.model.HistoryOutputRequest(
                    name="HISTORY_OUTPUT_" + step, createStepName=step,
                    variables=history_outputs)
                    # timeInterval=EVERY_TIME_INCREMENT)

        # Create and run the job
        job_name = "JOB_" + self.model_name
        mdb.Job(name=job_name, model=self.model_name, **kwargs)
        mdb.jobs[job_name].submit()
        mdb.jobs[job_name].waitForCompletion()

    def open_odb_file(self):
        """
        Open the output database file containing the simulation results.

        Returns
        -------
        None.

        """

        odb_file = "JOB_" + self.model_name + ".odb"
        if not os.path.isfile(odb_file):
            raise FileNotFoundError(
                "The file {} does not exist".format(odb_file))

        if odb_file not in session.odbs.keys():
            session.openOdb(name=odb_file)
        self._odb = session.odbs[odb_file]

    def close_odb_file(self):
        """
        Close the output database file.

        Returns
        -------
        None.

        """

        if self._odb is not None:
            self._odb.close()
        self._odb = None

    def export_to_vtk(self, path, field_outputs=None, node_sets=None,
                      element_sets=None, radius=None, clear=False):
        """
        Export the simulation results to VTK file format.

        The Python package 'Paraqus' needs to be installed to have this
        function available:

            https://github.com/InstituteOfMechanics/Paraqus/tree/main

        Parameters
        ----------
        path : str
            Path to the folder in which the files will be stored.
        field_outputs : list[str], optional
            The field outputs that shall be exported. Invariants can be
            requested in the form of a list element of length 2, e.g.
            ("S", "MISES"). The default is None.
        node_sets : list[str], optional
            The node sets that shall be exported. The default is None.
        element_sets : list[str], optional
            The element sets that shall be exported. The default is None.
        radius : dict[str, float], optional
            A map from element set names to the strut radius. It overrides
            or supplements the radii recorded during section assignment.
            The default is None, i.e. only the recorded radii are used.
        clear : bool, optional
            If True, the output folder will be cleared if it already
            exists. The default is False.

        Returns
        -------
        None.

        """

        if not PARAQUS_INSTALLED:
            raise ImportError("Paraqus is needed for the VTK export.")

        # Normalize the optional arguments (avoid mutable defaults)
        field_outputs = list(field_outputs) if field_outputs else []
        node_sets = list(node_sets) if node_sets else []
        element_sets = list(element_sets) if element_sets else []

        # Strut radii for the radius cell field
        strut_radius = self._section_radius
        if radius is not None:
            strut_radius.update(radius)

        def add_radius_field(instance_model):
            """
            Add radius as a field to truss and beam element sets.

            Parameters
            ----------
            instance_model : ParaqusModel
                The model that is converted to VTK format.

            Returns
            -------
            None.

            """

            if len(strut_radius) == 0:
                return

            element_radius = {}
            for name, value in strut_radius.items():
                set_name = (name if name.startswith("ELSET_")
                            else "ELSET_" + name)
                group = self.model_name + "." + set_name
                for tag in instance_model.elements.groups[group]:
                    element_radius[tag] = value
            instance_model.add_field("RADIUS", list(element_radius.keys()),
                                     list(element_radius.values()),
                                     "ELEMENTS", "SCALAR")

        # The strut element sets must be exported for the radius field
        for set_name in strut_radius.keys():
            if set_name not in element_sets:
                element_sets.append(set_name)

        if not os.path.isdir(path):
            os.makedirs(path)

        if self.odb is None:
            self.open_odb_file()

        # Initialize the reader and request fields and sets
        odb_file = "JOB_" + self.model_name + ".odb"
        reader = OdbReader(odb_path=odb_file,
                           model_name=self.model_name,
                           instance_names=[self.model_name])
        for output in field_outputs:
            if type(output) is not str and len(output) == 2:
                reader.add_field_export_request(output[0], invariant=output[1])
            else:
                reader.add_field_export_request(output)
        for name in element_sets:
            set_name = name if name.startswith("ELSET_") else "ELSET_" + name
            reader.add_set_export_request(set_name,
                                          set_type="ELEMENTS",
                                          instance_name=self.model_name)
        for name in node_sets:
            set_name = name if name.startswith("NSET_") else "NSET_" + name
            reader.add_set_export_request(set_name,
                                          set_type="NODES",
                                          instance_name=self.model_name)

        # Run over steps and frames for the export
        vtu_writer = BinaryWriter(path, clear_output_dir=clear)
        with CollectionWriter(vtu_writer, self.model_name) as writer:
            for step_name in self.odb.steps.keys():
                if step_name == "Initial":
                    continue
                for frame in reader.get_frame_indices(step_name):
                    instance_model = list(reader.read_instances(
                        step_name=step_name, frame_index=frame))[0]
                    add_radius_field(instance_model)
                    writer.write(instance_model)

        self.close_odb_file()

    def get_number_of_frames(self, step_name, close=False):
        """
        Return the number of frames for a given step.

        Parameters
        ----------
        step_name : str
            The name of the step.
        close : bool, optional
            If True, the ODB file will be closed after getting the
            number of frames. The default is False.

        Returns
        -------
        n_frames : int
            The number of frames available in the given step.

        """

        self.open_odb_file()
        n_frames = len(self.odb.steps[step_name].frames)
        self.close_odb_file()
        return n_frames

    def average_stress(self, step_name, frame=-1, share_boundary_struts=True):
        r"""
        Compute the averaged stress.

        In case of linear kinematics the average Cauchy stress is
        computed. In case of nonlinear kinematics, the first Piola-
        Kirchhoff stress is computed.

        Parameters
        ----------
        step_name : str
            Name of the step.
        frame : int, optional
            The frame index. The default is -1.
        share_boundary_struts : bool, optional
            If True, weight boundary struts by their periodic sharing
            fraction. The default is True.

        Returns
        -------
        stress : numpy.ndarray[float]
            The averaged stress tensor in 2d or 3d symbolic tensor
            notation.

        """

        self.open_odb_file()

        dim = 3 if self.model_dimension == THREE_D else 2

        # Compute the RVE volume
        bounding_box = np.array(self.get_bounding_box()).reshape(2, -1)
        delta = bounding_box[1] - bounding_box[0]
        if dim == 2:
            delta[-1] = 1.0
        volume = float(np.prod(delta))

        frame = self.odb.steps[step_name].frames[frame]

        # NFORC is stored as scalar components
        field_outputs = frame.fieldOutputs
        if "NFORC1" not in field_outputs.keys():
            raise RuntimeError("Field output 'NFORC' is not available.")

        odb_instance = self.odb.rootAssembly.instances[self.model_name]
        node_coords = {n.label: np.array(n.coordinates[0:dim])
                       for n in odb_instance.nodes}
        weights = (self._element_sharing_weights() if share_boundary_struts
                   else {})

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

        # Virial sum with the per-element sharing weights
        stress = np.zeros((dim, dim), dtype=float)
        for (element_label, node_label), force in forces.items():
            weight = weights.get(element_label, 1.0)
            stress += weight*np.outer(force, node_coords[node_label])
        stress *= -1.0/volume

        self.close_odb_file()
        return stress

    def _element_sharing_weights(self):
        r"""
        Compute the periodic sharing weight of every element.

        Returns
        -------
        weights : dict[int, float]
            A map from element label to sharing weight.

        """

        x_min, y_min, z_min, x_max, y_max, z_max = self._bounding_box
        planes = [(0, x_min), (0, x_max), (1, y_min), (1, y_max)]
        if self.model_dimension == THREE_D:
            planes += [(2, z_min), (2, z_max)]
        delta = max(x_max - x_min, y_max - y_min, z_max - z_min)
        tol = 1e-6*delta + 1e-12

        nodes = self.part.nodes
        weights = {}
        for element in self.part.elements:
            coords = []
            for index in element.connectivity:
                coords.append(np.array(nodes[index].coordinates))

            n_faces = 0
            for axis, value in planes:
                on_face = True
                for coord in coords:
                    if abs(coord[axis] - value) >= tol:
                        on_face = False
                        break
                if on_face:
                    n_faces += 1
            weights[element.label] = 0.5**n_faces
        return weights


class LoadHistory(object):
    """
    A time-dependent stress or strain load history for a homogenization.

    Parameters
    ----------
    time : ArrayLike
        The time data points of the load history
    history : Callable
        The stress or strain history function having a point of time as
        input and returns the load in full 3d tensor notation, i.e. in
        form of a (3, 3) array.

    Attributes
    ----------
    time_data : ndarray
        The time data points of the load history
    load_data : ndarray
        The load data points of the load history, i.e. the history
        evaluated at the time data points.

    Example
    -------
    >>> from abqhom.utils import voigt_to_tensor
    >>> from abqhom.abqmodel import LoadHistory
    >>>
    >>> # Strain history for a plane strain problem
    >>> def strain_path(time):
    >>>     strain = np.zeros(6)
    >>>     strain[0] = np.sin(np.pi*time)
    >>>     strain[1] = time**2
    >>>     strain[5] = np.sin(2*np.pi*time)
    >>>     strain = voigt_to_tensor(strain, "STRAIN")
    >>>     return strain
    >>>
    >>> time = np.linspace(0.0, 1.0, 50)
    >>> history = LoadHistory(time, strain_path)
    >>> initial_strain = history(0.0)

    """

    def __init__(
            self, time, history):
        self._time = np.array(time)
        self._history = history
        self._index = 0

        # Test return value of history
        self[self._index]

    @property
    def time_data(self):
        """
        The time data points of the load history

        Returns
        -------
        ndarray
            The time data points.

        """
        return self._time

    @property
    def load_data(self):
        """
        The load data points (stress or strain) of the load history.

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
        time, load = self.__next__()
        return time, load

    def __getitem__(self, i):
        time = self.time_data[i]
        load = np.array(self(time))

        if load.shape != (3, 3):
            raise ValueError(
                "The history function must return the load in full 3d tensor "
                "notation as a (3, 3) array.")

        return time, load

    def __call__(self, time):
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
        data : tuple[tuple[float, ndarray]]
            The (time, load) pairs.

        """

        data = []
        for t, v in self:
            data.append((t, np.array(v)))
        data = tuple(data)
        return data
