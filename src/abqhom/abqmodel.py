"""
Module containing the AbaqusRVE class that provides useful methods for
computational homogenization with Abaqus.

"""

import os
import numpy as np
from abqhom.utils import voigt_to_tensor, tensor_to_voigt
from abaqus import Mdb, mdb, session
from regionToolset import Region
from mesh import MeshNodeArray, MeshElementArray
from abaqusConstants import (THREE_D, TWO_D_PLANAR, DEFORMABLE_BODY, ON, UNSET,
                             CARTESIAN, N1_COSINES, DURING_ANALYSIS)
from abqhom.constants import (POINT_GROUPS, EDGE_GROUPS, FACE_GROUPS,
                              PERIODIC_PAIRS, ELEMENTS_ROTATIONAL_DOFS,
                              PARAQUS_INSTALLED)

if PARAQUS_INSTALLED:
    from paraqus.abaqus import OdbReader
    from paraqus.writers import BinaryWriter, CollectionWriter


class AbaqusRVE():
    """
    An RVE representation to set up a homogenization problem in Abaqus.

    Parameters
    ----------
    gmsh_RVE : GmshRVE
        The GmshRVE containing the mesh information of the RVE.
    input_file : str
        Path to the Abaqus input file containing the mesh information.

    Attributes
    ----------
    model_name : str
        The name of the model.
    model_dimension : abaqus.abaqusConstant
        The model dimension. The value can be THREE_D or TWO_D_PLANAR.
    model : abaqus.Model
        The internal Abaqus representation of the model.
        element sets
    part : abaqus.Part
        The internal Abaqus representation of the RVE part.
    instance : abaqus.PartInstance
        The internal representation of the RVE part in the assembly.
    odb : abaqus.Odb
        The Abaqus output database containing the analysis results.


    Example
    -------
    >>> # Assume there exists a GmshRVE and an input file
    >>> # Create an AbaqusRVE from that
    >>> unit_cell = AbaqusRVE(gmsh_RVE, input_file)
    >>>
    >>> # Assign material and sections
    >>> material = unit_cell.model.Material(name="Steel")
    >>> material.Elastic(table=[[210.0, 0.3]])
    >>> unit_cell.assign_solid_section("Steel", "Steel", "ELSET_ALL_ELEMENTS")
    >>>
    >>> # Generate assembly and step
    >>> unit_cell.model.StaticStep(name="Load", previous="Initial")
    >>>
    >>> # Apply periodic boundary conditions and run analysis
    >>> strain = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]  # Shear strain in xy-plane
    >>> unit_cell.apply_periodic_BC("Load", strain)
    >>> unit_cell.run_analysis(field_outputs=["NFORC", "U", "S", "E"])
    >>>
    >>> # Export the results to VTK file format
    >>> unit_cell.export_to_vtk(vtk_path, field_outputs=["U", "S", "E"])

    """

    def __init__(self, gmsh_RVE, input_file):
        if gmsh_RVE.model_dimension not in {2, 3}:
            raise RuntimeError("Only 2d and 3d models are supported.")
        dim = THREE_D if gmsh_RVE.model_dimension == 3 else TWO_D_PLANAR

        # Load the RVE into Abaqus
        Mdb()
        default_model = mdb.models.keys()[0]
        mdb.ModelFromInputFile(name=gmsh_RVE.model_name,
                               inputFileName=input_file)
        del mdb.models[default_model]

        self._model_name = gmsh_RVE.model_name
        self._model_dimension = dim
        self._model = mdb.models[self._model_name]
        self._part = self._model.parts[self._model_name]
        self._instance = None
        self._odb = None

        # Set the correct modeling space
        self._part.setValues(space=dim)

        # Add node sets
        if len(gmsh_RVE.node_sets) == 0:
            gmsh_RVE.build_node_sets()
        self._add_node_sets(gmsh_RVE.node_sets)

        # Add element sets
        if len(gmsh_RVE.element_sets) == 0:
            gmsh_RVE.build_element_sets()
        self._add_element_sets(gmsh_RVE.element_sets)

        # Add part instance to the assembly
        assembly = self.model.rootAssembly
        assembly.DatumCsysByDefault(CARTESIAN)
        self._instance = assembly.Instance(name=self._model_name,
                                           part=self._part,
                                           dependent=ON)

    @property
    def model_name(self):
        """The name of the Abaqus model."""
        return self._model_name

    @property
    def model_dimension(self):
        """The geometric dimensionality of the Abaqus model."""
        return self._model_dimension

    @property
    def model(self):
        """A shortcut to the RVE Abaqus model."""
        return self._model

    @property
    def part(self):
        """A shortcut to the RVE part in the Abaqus model."""
        return self._part

    @property
    def instance(self):
        """A shortcut to the RVE part instance in the Abaqus model."""
        return self._instance

    @property
    def odb(self):
        """A shortcut to the output data base file resulting from the model."""
        return self._odb

    def _add_node_sets(self, node_sets):
        """Add all relevant node sets to the Abaqus model."""

        for set_name, set_items in node_sets.items():
            set_name = "NSET_" + str(set_name)
            set_items = [int(i) for i in set_items]

            # Sets must be unsorted to not break ordering for periodic BCs
            if len(set_items) > 0:
                self.part.SetFromNodeLabels(name=set_name,
                                            nodeLabels=set_items,
                                            unsorted=True)

    def _add_element_sets(self, element_sets):
        """Add all relevant element sets to the Abaqus model."""

        for set_name, set_items in element_sets.items():
            set_name = "ELSET_" + str(set_name)
            set_items = [int(i) for i in set_items]

            if len(set_items) > 0:
                self.part.SetFromElementLabels(name=set_name,
                                               elementLabels=set_items)

    def _create_reference_point(self, direction):
        """Create a reference points to apply the BCs in one direction."""

        coords = [0.0, 0.0, 0.0]
        name = "REFERENCE_POINT_X{}".format(direction)
        part = self.model.Part(dimensionality=self.model_dimension,
                               name=name, type=DEFORMABLE_BODY)
        point = self.model.parts[name].ReferencePoint(point=coords)
        self.model.parts[name].Set(
            name="POINT_BOUNDARY_CONDITIONS",
            referencePoints=[self.model.parts[name]
                             .referencePoints[point.id]])
        self.model.rootAssembly.Instance(name=name, part=part,
                                         dependent=ON)

    def _get_bounding_box(self):
        """Compute the bounding box of the Abaqus model."""

        if self.model_dimension == THREE_D:
            x_min, y_min, z_min = (self.instance.sets["NSET_POINT_1"].nodes[0]
                                   .coordinates)
            x_max, y_max, z_max = (self.instance.sets["NSET_POINT_7"].nodes[0]
                                   .coordinates)
        else:  # TWO_D_PLANAR
            x_min, y_min, z_min = (self.instance.sets["NSET_POINT_1"].nodes[0]
                                   .coordinates)
            x_max, y_max, z_max = (self.instance.sets["NSET_POINT_3"].nodes[0]
                                   .coordinates)

        return x_min, y_min, z_min, x_max, y_max, z_max

    def _collect_element_types(self):
        """Create a set containing the names of all occuring element types."""
        all_elements = self.part.sets["ELSET_ALL_ELEMENTS"].elements
        element_types = {str(element.type) for element in all_elements}
        return element_types

    def apply_uniform_traction_BC(self, step_name, stress):
        """
        Apply uniform traction boundary conditions to the RVE model.

        Parameters
        ----------
        step_name : str
            The name of the step in which the boundary conditions are applied.
        stress : numpy.ndarray[float]
            The macroscopic stress tensor in full 3d Voigt notation, i.e. the
            order of the components is [xx, yy, zz, yz, xz, xy].

        Returns
        -------
        None.

        """

        raise NotImplementedError(
            "Uniform traction boundary conditions are not supported yet.")

    def apply_affine_displacement_BC(self, step_name, strain):
        """
        Apply affine displacement boundary conditions to the RVE model.

        Parameters
        ----------
        step_name : str
            The name of the step in which the boundary conditions are applied.
        strain : numpy.ndarray[float]
            The macroscopic strain tensor in full 3d Voigt notation, i.e. the
            order of the components is [xx, yy, zz, 2yz, 2xz, 2xy].

        Returns
        -------
        None.

        """

        if self._collect_element_types() <= ELEMENTS_ROTATIONAL_DOFS:
            raise NotImplementedError(
                "Affine displacement boundary conditions are not supported "
                + "yet for elements with rotational degrees of freedom. "
                + "Please use periodic boundary conditions.")

        if not np.array(strain).shape == (6,):
            raise ValueError(
                "Strain must be provided in full 3d Voigt notation.")
        strain = voigt_to_tensor(np.array(strain), "STRAIN")

        # Apply Dirichlet boundary conditions on all boundary nodes
        boundary_nodes = self.instance.sets["NSET_BOUNDARY_NODES"]
        for node in boundary_nodes.nodes:
            region = Region(nodes=MeshNodeArray(nodes=[node]))
            x = np.array(node.coordinates)

            u = np.dot(strain, x)
            u1, u2 = u[0], u[1]
            u3 = u[2] if self.model_dimension == THREE_D else UNSET

            name = "DIRICHLET_{}".format(node.label)
            if name not in self.model.boundaryConditions:
                self.model.DisplacementBC(name=name,
                                          createStepName=step_name,
                                          region=region,
                                          u1=u1, u2=u2, u3=u3)
            else:
                bc = self.model.boundaryConditions[name]
                bc.setValuesInStep(stepName=step_name, u1=u1, u2=u2, u3=u3)

    def apply_periodic_BC(self, step_name, strain):
        """
        Apply periodic boundary conditions to the RVE model.

        Parameters
        ----------
        step_name : str
            The name of the step in which the boundary conditions are applied.
        strain : numpy.ndarray[float]
            The macroscopic strain tensor in full 3d Voigt notation, i.e. the
            order of the components is [xx, yy, zz, 2yz, 2xz, 2xy].

        Returns
        -------
        None.

        """

        if not np.array(strain).shape == (6,):
            raise ValueError(
                "Strain must be provided in full 3d Voigt notation.")
        strain = voigt_to_tensor(np.array(strain), "STRAIN")

        def apply_dirichlet_BC_to_reference_point(direction):
            """Apply BCs to the reference points for periodic connections."""

            # Compute RVE size
            bounding_box = np.array(self._get_bounding_box())
            bounding_box = bounding_box.reshape(2, -1)
            deltas = bounding_box[1] - bounding_box[0]

            # Apply the boundary conditions
            assembly = self.model.rootAssembly
            point_name = "REFERENCE_POINT_X{}".format(direction)
            set_name = "POINT_BOUNDARY_CONDITIONS"
            region = assembly.instances[point_name].sets[set_name]
            x = np.array([0.0, 0.0, 0.0])
            x[direction - 1] = deltas[direction - 1]
            u = np.dot(strain, x)
            u1, u2 = u[0], u[1]
            u3 = u[2] if self.model_dimension == THREE_D else UNSET
            self.model.DisplacementBC(name="DIRICHLET_" + point_name,
                                      createStepName=step_name,
                                      region=region,
                                      u1=u1, u2=u2, u3=u3)

        def apply_periodic_rotations():
            """Apply periodic boundary conditions to rotational DOFs."""

            # Make a list of all periodically connected pairs of groups
            all_pairs = {i: PERIODIC_PAIRS[i] for i in PERIODIC_PAIRS}
            for i in range(len(POINT_GROUPS) - 1):
                all_pairs[POINT_GROUPS[i]] = POINT_GROUPS[i + 1]

            # Apply equation constraints to rotational degrees of freedom
            for master, slave in all_pairs.items():
                if ("NSET_" + master not in self.instance.sets
                        or "NSET_" + slave not in self.instance.sets):
                    continue

                dofs = [4, 5, 6] if self.model_dimension == THREE_D else [6]
                for dof in dofs:
                    base_name = self.model_name + ".NSET_"
                    terms = [[1.0, base_name + slave, dof],
                             [-1.0, base_name + master, dof]]

                    self.model.Equation(name="PBC_{}_{}_{}".format(
                        master, slave, dof), terms=terms)

        # Define some shortcuts
        RVE_sets = self.instance.sets
        dim = 3 if self.model_dimension == THREE_D else 2

        # Apply Dirichlet BC to corner points
        for point in POINT_GROUPS:
            if "NSET_" + point not in RVE_sets:
                continue
            region = RVE_sets["NSET_" + point]
            x = np.array(region.nodes[0].coordinates)

            u = np.dot(strain, x)
            u1, u2 = u[0], u[1]
            u3 = u[2] if dim == 3 else UNSET
            self.model.DisplacementBC(name="DIRICHLET_" + point,
                                      createStepName=step_name,
                                      region=region,
                                      u1=u1, u2=u2, u3=u3)

        # Add periodic coupling to edges and faces
        for master, slave in PERIODIC_PAIRS.items():
            if ("NSET_" + master not in RVE_sets
                    or "NSET_" + slave not in RVE_sets):
                continue

            x_minus = RVE_sets["NSET_" + master].nodes[0].coordinates[0:dim]
            x_plus = RVE_sets["NSET_" + slave].nodes[0].coordinates[0:dim]
            distance = np.array(x_plus) - np.array(x_minus)
            direction = np.argmax(np.abs(distance))

            # Create a reference point for the RHS of the equation
            reference_point = "REFERENCE_POINT_X{}".format(direction + 1)
            if reference_point not in self.model.parts:
                self._create_reference_point(direction + 1)
                apply_dirichlet_BC_to_reference_point(direction + 1)

            for dof in range(1, dim + 1):
                base_name = self.model_name + ".NSET_"
                ref_set_name = reference_point + ".POINT_BOUNDARY_CONDITIONS"
                terms = [[1.0, base_name + slave, dof],
                         [-1.0, base_name + master, dof],
                         [-1.0*np.sign(distance[direction]),
                          ref_set_name, dof]]

                self.model.Equation(name="PBC_{}_{}_{}".format(
                    master, slave, dof), terms=terms)

        # Apply constraints also on rotational degrees of freedom
        if self._collect_element_types() <= ELEMENTS_ROTATIONAL_DOFS:
            apply_periodic_rotations()

    def run_analysis(self, field_outputs=["NFORC"], history_outputs=[],
                     **kwargs):
        """
        Create the job and run the finite element analysis in Abaqus.

        Parameters
        ----------
        field_outputs : list[str], optional
            The field outputs that shall be computed. The default is ["NFORC"].
        history_outputs : list[str], optional
            The history outputs that shall be computed. The default is [].
        **kwargs :
            Additional arguments that will be passed to the Abaqus job. See
            the documentation of the Abaqus Job class for details.

        Returns
        -------
        None.

        """
        # Remove defaults
        if "F-Output-1" in self.model.fieldOutputRequests:
            del self.model.fieldOutputRequests["F-Output-1"]
        if "H-Output-1" in self.model.historyOutputRequests:
            del self.model.historyOutputRequests["H-Output-1"]

        # Request the field outputs
        if len(field_outputs) > 0:
            for step in self.model.steps.keys():
                if step == "Initial":
                    continue
                self.model.FieldOutputRequest(name="FIELD_OUTPUT_" + step,
                                              createStepName=step,
                                              variables=field_outputs)

        # Request the history outputs
        if len(history_outputs) > 0:
            for step in self.model.steps.keys():
                if step == "Initial":
                    continue
                self.model.HistoryOutputRequest(name="HISTORY_OUTPUT_" + step,
                                                createStepName=step,
                                                variables=field_outputs)

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

    def assign_solid_section(self, section_name, material_name,
                             element_set_name):
        """
        Assign a solid material section to a given element set.

        Parameters
        ----------
        section_name : str
            The name of the section to be created.
        material_name : str
            The name of the material that is associated with the section. The
            material should exist already.
        element_set_name : str
            The name of the element set that is assigned the section.

        Returns
        -------
        None.

        """

        self.model.HomogeneousSolidSection(name=section_name,
                                           material=material_name)
        self._assign_general_section(section_name, element_set_name)

    def assign_beam_section(self, section_name, material_name,
                            element_set_name, radius):
        """
        Assign a beam section with circular profile to a given element set.

        The orientation of the beams will be created automatically.

        Parameters
        ----------
        section_name : str
            The name of the section to be created.
        material_name : str
            The name of the material that is associated with the section. The
            material should exist already.
        element_set_name : str
            The name of the element set that is assigned the section.
        radius : float
            The radius of the circular cross-sectional area.

        Returns
        -------
        None.

        """

        self.model.CircularProfile(name=section_name, r=radius)
        self.model.BeamSection(name=section_name, integration=DURING_ANALYSIS,
                               profile=section_name, material=material_name)
        self._assign_general_section(section_name, element_set_name)
        self._assign_beam_section_orientation(element_set_name)

    def assign_shell_section(self, section_name, material_name,
                             element_set_name):
        """
        This method is not supported yet.

        Parameters
        ----------
        section_name : str
            The name of the section to be created.
        material_name : str
            The name of the material that is associated with the section. The
            material should exist already.
        element_set_name : str
            The name of the element set that is assigned the section.

        Returns
        -------
        None.

        """
        raise NotImplementedError("Shell sections are not supported yet.")

    def assign_truss_section(self, section_name, material_name,
                             element_set_name, area):
        """
        Assign a truss section to a given element set.

        Parameters
        ----------
        section_name : str
            The name of the section to be created.
        material_name : str
            The name of the material that is associated with the section. The
            material should exist already.
        element_set_name : str
            The name of the element set that is assigned the section.
        area : float
            The cross-sectional area of the beam.

        Returns
        -------
        None.

        """

        self.model.TrussSection(name=section_name, material=material_name,
                                area=area)
        self._assign_general_section(section_name, element_set_name)

    def _assign_general_section(self, section_name, element_set_name):
        """Assign a material section to a given element set."""
        if element_set_name not in self.part.sets:
            element_set_name = "ELSET_{}".format(element_set_name)
        region = self.part.sets[element_set_name]
        self.part.SectionAssignment(region=region, sectionName=section_name)

    def _assign_beam_section_orientation(self, element_set_name):
        """Compute the orientation of beam elements with circular profile."""
        if element_set_name not in self.part.sets:
            element_set_name = "ELSET_{}".format(element_set_name)

        for element in self.part.sets[element_set_name].elements:
            # Compute unit tangent vector, i.e. the beam axis direction
            conn = element.connectivity
            node_start = self.part.nodes[conn[0]]
            node_end = self.part.nodes[conn[-1]]
            coords_start = np.array(node_start.coordinates)
            coords_end = np.array(node_end.coordinates)
            tangent = coords_end - coords_start
            tangent /= np.linalg.norm(tangent)

            # Compute a vector orthogonal to the beam axis
            free_index = 0
            all_indices = [0, 1, 2]
            while tangent[free_index] == 0:
                free_index += 1
            all_indices.remove(free_index)
            normal = np.ones_like(tangent)
            normal[free_index] = -1/tangent[free_index]*(
                tangent[all_indices[0]] + tangent[all_indices[1]])
            normal /= np.linalg.norm(normal)

            assert np.dot(tangent, normal) == 0

            # Assign orientation to the element
            region = Region(elements=MeshElementArray(elements=[element]))
            self.part.assignBeamSectionOrientation(region=region,
                                                   method=N1_COSINES,
                                                   n1=normal)

    def export_to_vtk(self, path, field_outputs=[], sets=[], radius={},
                      clear=False):
        """
        Export the simulation results to VTK file format.

        The Python package 'Paraqus' needs to be installed to have this
        function available.

        Parameters
        ----------
        path : str
            Path to the folder in which the files will be stored.
        field_outputs : list[str], optional
            The field outputs that shall be exported. The default is [].
        sets : list[str], optional
            The sets that shall be exported. Invariants can be requested
            in the form of a list element of length 2, e.g. ("S", "MISES").
            The default is [].
        diameters : dict{str, float}
            DESCRIPTION
        clear : bool, optional
            If True, the output folder will be cleared if it already exists.
            The default is False.

        Returns
        -------
        None.

        """
        if not PARAQUS_INSTALLED:
            raise ImportError("Paraqus is needed for the VTK export.")

        if not os.path.isdir(path):
            os.makedirs(path)

        if self.odb is None:
            self.open_odb_file()

        # Initialize the reader and request fields and sets
        odb_file = "JOB_" + self.model_name + ".odb"
        reader = OdbReader(odb_path=odb_file,
                           model_name=self.model_name,
                           instance_names=[self.model_name])
        for i in field_outputs:
            if len(i) > 1:
                reader.add_field_export_request(i[0], invariant=i[1])
            else:
                reader.add_field_export_request(i)

        for i in radius.keys():
            if i not in sets:
                sets.append(i)
        for i in sets:
            if i.startswith("NSET"):
                position = "NODES"
            elif i.startswith("ELSET"):
                position = "ELEMENTS"
            reader.add_set_export_request(i,
                                          set_type=position,
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

                    # Add diameter fields if assigned
                    element_tags = []
                    radius_values = []
                    for i, j in radius.items():
                        group_elements = instance_model.elements.groups[
                            self.model_name + "." + i]
                        element_tags.extend(group_elements)
                        radius_values.extend(
                            np.full(group_elements.shape, j))
                    instance_model.add_field("RADIUS", element_tags,
                                             radius_values, "ELEMENTS",
                                             "SCALAR")

                    writer.write(instance_model)

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
        self.open_odb_file()
        n_frames = len(self.odb.steps[step_name].frames)
        return n_frames

    def compute_averaged_stress(self, step_name, frame=-1):
        """
        Compute the averaged stress.

        NOTE: This function should be re-implemented as soon as we know how
        to efficiently compute averages for elasto-plastic fields such as
        the plastic strain.

        Parameters
        ----------
        step_name : str
            Name of the step.
        frame : int, optional
            The frame index. The default is -1.

        Returns
        -------
        stress : numpy.ndarray[flot]
            The averaged stress tensor in 2d or 3d Voigt notation.

        """
        self.open_odb_file()

        dim = 3 if self.model_dimension == THREE_D else 2
        bounding_box = np.array(self._get_bounding_box())
        bounding_box = bounding_box.reshape(2, -1)
        deltas = bounding_box[1] - bounding_box[0]
        deltas = deltas[0:dim]

        # Get internal forces
        frame = self.odb.steps[step_name].frames[frame]
        f_int_1 = frame.fieldOutputs["NFORC1"]
        f_int_2 = frame.fieldOutputs["NFORC2"]
        f_int_3 = (frame.fieldOutputs["NFORC3"]
                   if self.model_dimension == THREE_D else f_int_1)

        # Initialize averaged stress tensor
        stress = np.zeros((dim, dim), dtype=float)

        # Loop over all boundary nodes
        boundary_sets = POINT_GROUPS + EDGE_GROUPS
        if dim == 3:
            boundary_sets += FACE_GROUPS
        odb_instance = self.odb.rootAssembly.instances[self.model_name]
        for set_name in boundary_sets:
            if "NSET_" + set_name not in odb_instance.nodeSets:
                continue
            region = odb_instance.nodeSets["NSET_" + set_name]
            f_int_1_subset = f_int_1.getSubset(region=region).values
            f_int_2_subset = f_int_2.getSubset(region=region).values
            f_int_3_subset = f_int_3.getSubset(region=region).values

            # Collect force values for each node
            f_int_1_map, f_int_2_map, f_int_3_map = {}, {}, {}
            for i, j, k in zip(f_int_1_subset, f_int_2_subset, f_int_3_subset):

                f_int_1_map[i.nodeLabel] = (f_int_1_map.get(i.nodeLabel,
                                                            0.0) - i.data)
                f_int_2_map[j.nodeLabel] = f_int_2_map.get(j.nodeLabel,
                                                           0.0) - j.data
                f_int_3_map[k.nodeLabel] = f_int_3_map.get(k.nodeLabel,
                                                           0.0) - k.data

            # Update averaged stress tensor
            for node in region.nodes:
                x = node.coordinates[0:dim]
                f_int = np.array([f_int_1_map[node.label],
                                  f_int_2_map[node.label],
                                  f_int_3_map[node.label]])
                f_int = f_int[0:dim]
                stress += np.outer(f_int, x)

        # Divide the stress by the RVE volume
        stress /= np.prod(deltas)
        stress = 0.5*(stress + stress.T)
        stress = tensor_to_voigt(stress, "STRESS")

        self.close_odb_file()
        return stress
