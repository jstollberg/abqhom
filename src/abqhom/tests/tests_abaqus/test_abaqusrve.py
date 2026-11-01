# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-only

"""
Granular tests for AbaqusRVE that exercise individual pieces of the API
without running a finite element analysis: construction, set transfer,
the removal of redundant periodic elements, the rigid-body anchor
selection and the recorded section radii. These need the Abaqus kernel.

"""
import unittest

import numpy as np

import abqhom.abqmodel as abqmodel
from abqhom.utils import voigt_to_tensor
from abqhom.gmshmodel import GmshRVE
from abqhom.abqmodel import AbaqusRVE, LoadHistory
from abqhom.tests import fixtures
from abqhom.tests.base import AbaqusTestCase

from abaqusConstants import (THREE_D, TWO_D_PLANAR, ON, OFF, DURING_ANALYSIS,
                             BEFORE_ANALYSIS, TIME, FREQUENCY)


def _uniaxial_history():
    """A trivial uniaxial (symmetric) strain history for BC tests."""
    def strain_history(time):
        strain = np.zeros(6)
        strain[0] = 1.0e-3
        return time*voigt_to_tensor(strain, "STRAIN")
    return LoadHistory(np.array([0.0, 1.0]), strain_history)


def _build(builder, element_size=1.0, custom_types=None):
    """Build, mesh, and transfer a fixture to an AbaqusRVE."""
    name = builder()
    rve = GmshRVE(name)
    rve.add_periodic_mesh_constraints()
    rve.generate_mesh(element_size=element_size)
    return AbaqusRVE(rve, custom_types or {})


class TestLoadHistory(unittest.TestCase):

    @staticmethod
    def _history(time):
        """Return a simple tensor-valued history for validation tests."""
        return np.eye(3)*time

    def test_rejects_empty_time_data(self):
        """A load history requires at least one time value."""
        self.assertRaises(ValueError, LoadHistory, [], self._history)

    def test_rejects_multidimensional_time_data(self):
        """Time data must be a one-dimensional sequence."""
        self.assertRaises(ValueError, LoadHistory, [[0.0, 1.0]],
                          self._history)

    def test_rejects_nonfinite_time_data(self):
        """Time data may not contain nonfinite values."""
        self.assertRaises(ValueError, LoadHistory, [0.0, np.inf],
                          self._history)

    def test_rejects_nonincreasing_time_data(self):
        """Time values must be strictly increasing."""
        self.assertRaises(ValueError, LoadHistory, [0.0, 0.0], self._history)
        self.assertRaises(ValueError, LoadHistory, [1.0, 0.0], self._history)

    def test_rejects_noncallable_history(self):
        """The prescribed history must be callable."""
        self.assertRaises(TypeError, LoadHistory, [0.0, 1.0], None)

    def test_rejects_history_with_invalid_tensor_shape(self):
        """The history function must return a full second-order tensor."""
        self.assertRaises(ValueError, LoadHistory, [0.0, 1.0],
                          lambda time: np.zeros(6))

    def test_iteration_properties_and_table_are_consistent(self):
        """All public views expose the same evaluated load history."""
        times = np.array([0.0, 0.5, 1.0])
        history = LoadHistory(times, self._history)

        self.assertTrue(np.array_equal(history.time_data, times))
        expected = np.array([self._history(time) for time in times])
        self.assertTrue(np.allclose(history.load_data, expected))

        table = history.as_table()
        self.assertEqual(len(table), len(times))
        for index, (time, load) in enumerate(table):
            self.assertEqual(time, times[index])
            self.assertTrue(np.allclose(load, expected[index]))


class TestConstruction(AbaqusTestCase):

    def test_vtk_reader_uses_current_job_name(self):
        """VTK export reads the most recently submitted repeated job."""
        class Reader(object):
            def __init__(self, odb_path, **kwargs):
                self.odb_path = odb_path

        original_reader = getattr(abqmodel, "OdbReader", None)
        abqmodel.OdbReader = Reader
        try:
            unit = object.__new__(AbaqusRVE)
            unit._model_name = "REPEATED"
            unit._job_name = "JOB_REPEATED_6"
            reader = unit._create_vtk_reader([], [], [])
            self.assertEqual(reader.odb_path, "JOB_REPEATED_6.odb")
        finally:
            if original_reader is None:
                delattr(abqmodel, "OdbReader")
            else:
                abqmodel.OdbReader = original_reader

    def test_missing_odb_raises_io_error(self):
        """Opening a nonexistent output database raises an I/O error."""
        unit = object.__new__(AbaqusRVE)
        unit._model_name = "ODB_THAT_DOES_NOT_EXIST"
        unit._odb = None
        self.assertRaises(IOError, unit.open_odb_file)

    def test_model_dimension_3d(self):
        """A spatial Gmsh cell becomes a three-dimensional Abaqus RVE."""
        unit = _build(lambda: fixtures.simple_cubic("A", 2.0))
        self.assertEqual(unit.model_dimension, THREE_D)
        self.assertIsNotNone(unit.part)
        self.assertIsNotNone(unit.instance)

    def test_model_dimension_2d(self):
        """A planar Gmsh cell becomes a planar Abaqus RVE."""
        unit = _build(lambda: fixtures.square_lattice("B", 2.0))
        self.assertEqual(unit.model_dimension, TWO_D_PLANAR)

    def test_bounding_box_matches_geometry(self):
        """The transferred Abaqus RVE retains the Gmsh bounding box."""
        unit = _build(lambda: fixtures.simple_cubic("C", 2.0))
        self.assertTrue(np.allclose(unit.get_bounding_box(),
                                    [-1, -1, -1, 1, 1, 1]))

    def test_volume_matches_3d_geometry(self):
        """RVE volume equals the cuboidal cell volume in three dimensions."""
        unit = _build(lambda: fixtures.simple_cubic("C3", 2.0))
        self.assertAlmostEqual(unit.get_volume(), 8.0)

    def test_volume_is_area_for_2d_geometry(self):
        """RVE volume denotes the cell area for a planar model."""
        unit = _build(lambda: fixtures.square_lattice("C2", 2.0))
        self.assertAlmostEqual(unit.get_volume(), 4.0)

    def test_frame_time_includes_time_accumulated_before_step(self):
        """Frame time includes the total time accumulated before its step."""
        class Frame(object):
            def __init__(self, value):
                self.frameValue = value

        class Step(object):
            domain = TIME
            totalTime = 2.0
            frames = [Frame(0.0), Frame(0.75)]

        class Odb(object):
            steps = {"LOAD": Step()}

        unit = object.__new__(AbaqusRVE)
        unit._odb = None
        unit.open_odb_file = lambda: setattr(unit, "_odb", Odb())
        unit.close_odb_file = lambda: setattr(unit, "_odb", None)

        self.assertAlmostEqual(unit.get_frame_time("LOAD", -1), 2.75)
        self.assertIsNone(unit.odb)

    def test_frame_time_rejects_non_time_domain_step(self):
        """Frequency-domain frames do not define a total analysis time."""
        class Step(object):
            domain = FREQUENCY

        class Odb(object):
            steps = {"FREQUENCY": Step()}

        unit = object.__new__(AbaqusRVE)
        unit._odb = None
        unit.open_odb_file = lambda: setattr(unit, "_odb", Odb())
        unit.close_odb_file = lambda: setattr(unit, "_odb", None)

        self.assertRaises(ValueError, unit.get_frame_time, "FREQUENCY")
        self.assertIsNone(unit.odb)

    def test_number_of_frames_closes_odb(self):
        """A method-owned ODB is closed after reading the frame count."""
        class Step(object):
            frames = [object(), object()]

        class Odb(object):
            steps = {"LOAD": Step()}

        unit = object.__new__(AbaqusRVE)
        unit._odb = None
        unit.open_odb_file = lambda: setattr(unit, "_odb", Odb())
        unit.close_odb_file = lambda: setattr(unit, "_odb", None)

        self.assertEqual(unit.get_number_of_frames("LOAD"), 2)
        self.assertIsNone(unit.odb)

    def test_average_stress_closes_odb_after_error(self):
        """A method-owned ODB is also closed along exceptional paths."""
        class Frame(object):
            fieldOutputs = {}

        class Step(object):
            frames = [Frame()]

        class Odb(object):
            steps = {"LOAD": Step()}

        unit = object.__new__(AbaqusRVE)
        unit._model_dimension = THREE_D
        unit._bounding_box = (0.0, 0.0, 0.0, 1.0, 1.0, 1.0)
        unit._odb = None
        unit.open_odb_file = lambda: setattr(unit, "_odb", Odb())
        unit.close_odb_file = lambda: setattr(unit, "_odb", None)

        self.assertRaises(RuntimeError, unit.average_stress, "LOAD")
        self.assertIsNone(unit.odb)

    def test_node_and_element_sets_are_transferred(self):
        """Physical and generated sets are transferred into Abaqus."""
        unit = _build(lambda: fixtures.simple_cubic("D", 2.0))
        self.assertIn("ELSET_STRUTS", unit.part.sets)
        self.assertIn("NSET_ALL_NODES", unit.part.sets)
        self.assertIn("NSET_POINT_1", unit.part.sets)

    def test_secondary_elements_are_removed_but_nodes_are_retained(self):
        """Reduction removes secondary elements without deleting nodes."""
        fixtures.simple_cubic("D2", length=2.0)
        rve = GmshRVE("D2")
        rve.add_periodic_mesh_constraints()
        rve.generate_mesh(element_size=1.0)
        rve.build_element_sets()
        all_nodes, _ = rve.get_nodes()
        remaining_elements = (len(rve.element_sets["ALL_ELEMENTS"])
                              - len(rve.element_sets["SECONDARY_ELEMENTS"]))

        unit = AbaqusRVE(rve)

        self.assertEqual(len(unit.part.nodes), len(all_nodes))
        self.assertEqual(len(unit.part.elements), remaining_elements)
        self.assertEqual(len(unit.part.sets["ELSET_STRUTS"].elements),
                         remaining_elements)

    def test_secondary_elements_can_be_retained(self):
        """The unreduced periodic representation remains available."""
        fixtures.simple_cubic("D3", length=2.0)
        rve = GmshRVE("D3")
        rve.add_periodic_mesh_constraints()
        rve.generate_mesh(element_size=1.0)
        rve.build_element_sets()

        unit = AbaqusRVE(rve, remove_secondary=False)

        self.assertEqual(len(unit.part.elements),
                         len(rve.element_sets["ALL_ELEMENTS"]))
        self.assertEqual(
            len(unit.part.sets["ELSET_SECONDARY_ELEMENTS"].elements),
            len(rve.element_sets["SECONDARY_ELEMENTS"]))

    def test_orphaned_structural_nodes_keep_periodic_rotations(self):
        """Removed elements do not erase their nodes' rotational traces."""
        fixtures.simple_cubic("D4", length=2.0)
        rve = GmshRVE("D4")
        rve.add_periodic_mesh_constraints()
        rve.generate_mesh(element_size=1.0)
        unit = AbaqusRVE(rve)

        retained = set()
        for element in unit.part.elements:
            for index in element.connectivity:
                retained.add(unit.part.nodes[index].label)
        orphaned = {node.label for node in unit.part.nodes} - retained

        self.assertTrue(orphaned)
        self.assertTrue(orphaned.issubset(
            set(unit._unreduced_node_element_types)))

        unit._periodic_boundary_rotations()
        constrained = set()
        for name in unit.part.sets.keys():
            if name.endswith("_DOF_4"):
                constrained.update(node.label
                                   for node in unit.part.sets[name].nodes)
        self.assertTrue(orphaned.issubset(constrained))


class TestRigidBodyAnchor(AbaqusTestCase):

    def test_anchor_is_a_corner_for_a_cubic_cell(self):
        """Automatic rigid-body anchoring selects a primary corner."""
        unit = _build(lambda: fixtures.simple_cubic("K", 2.0))
        region, node = unit.select_rigid_body_anchor()
        corner_labels = [n.label for n
                         in unit.instance.sets["NSET_POINT_1"].nodes]
        self.assertIn(node.label, corner_labels)

    def test_anchor_is_not_a_periodic_secondary_node(self):
        """The auto-selected anchor must be a primary node so it can
        carry a Dirichlet BC without over-constraining a periodic
        equation."""
        unit = _build(lambda: fixtures.plus_lattice_2d("L", 2.0))
        _, node = unit.select_rigid_body_anchor()
        from abqhom.constants import PERIODIC_PAIRS
        secondary_labels = set()
        for secondary in PERIODIC_PAIRS.keys():
            set_name = "NSET_" + secondary
            if set_name in unit.instance.sets:
                for n in unit.instance.sets[set_name].nodes:
                    secondary_labels.add(n.label)
        self.assertNotIn(node.label, secondary_labels)

    def test_user_defined_anchor_is_used(self):
        """A user-specified node set controls rigid-body anchoring."""

        # The BCC cell has an interior center node, so INNER_NODES is not
        # empty and can be used as a user-defined anchor.
        unit = _build(lambda: fixtures.bcc("M", 2.0))
        _, node = unit.select_rigid_body_anchor(fixed_node="INNER_NODES")
        inner_labels = [n.label for n
                        in unit.instance.sets["NSET_INNER_NODES"].nodes]
        self.assertIn(node.label, inner_labels)

    def test_anchor_is_not_an_embedded_node(self):
        """Automatic anchoring excludes translation-eliminated nodes."""
        unit = _build(lambda: fixtures.simple_cubic("M2", 2.0),
                      element_size=1.0)
        embedded = unit.instance.sets["NSET_POINT_1"].nodes[0].label
        unit._embedded_region_nodes = lambda: {embedded}
        _, node = unit.select_rigid_body_anchor()
        self.assertNotEqual(node.label, embedded)


class TestPeriodicBoundaryConditions(AbaqusTestCase):

    def test_one_sided_embedded_pair_is_rejected(self):
        """A periodic pair cannot contain only one eliminated translation."""
        unit = _build(lambda: fixtures.plus_lattice_2d("PBC_EMB", 2.0))
        unit.model.StaticStep(name="LOAD", previous="Initial")
        embedded = unit.part.sets["NSET_EDGE_2"].nodes[0].label
        unit._embedded_region_nodes = lambda: {embedded}
        self.assertRaises(RuntimeError, unit.apply_periodic_BC,
                          "LOAD", _uniaxial_history())


class TestSectionRadius(AbaqusTestCase):

    def test_beam_section_records_the_radius(self):
        """Beam-section assignment records its circular radius."""
        unit = _build(lambda: fixtures.simple_cubic("N", 2.0))
        material = unit.model.Material(name="MAT")
        material.Elastic(table=[[1.0, 0.3]])
        unit.assign_beam_section("STRUTS", "MAT", 0.07)
        self.assertAlmostEqual(unit._section_radius["STRUTS"], 0.07)

    def test_truss_section_records_the_equivalent_radius(self):
        """Truss assignment records the radius equivalent to its area."""
        unit = _build(lambda: fixtures.simple_cubic("O", 2.0))
        material = unit.model.Material(name="MAT")
        material.Elastic(table=[[1.0, 0.3]])
        area = 0.05
        unit.assign_truss_section("STRUTS", "MAT", area)
        self.assertAlmostEqual(unit._section_radius["STRUTS"],
                               (area/np.pi)**0.5)

    def test_beam_section_is_integrated_during_analysis(self):
        """By default the beam section is integrated during the analysis."""
        unit = _build(lambda: fixtures.simple_cubic("P", 2.0))
        material = unit.model.Material(name="MAT")
        material.Elastic(table=[[1.0, 0.3]])
        unit.assign_beam_section("STRUTS", "MAT", 0.07)
        self.assertEqual(unit.model.sections["STRUTS"].integration,
                         DURING_ANALYSIS)

    def test_beam_section_can_be_pre_integrated(self):
        """pre_integrate=True pre-integrates the beam section."""
        unit = _build(lambda: fixtures.simple_cubic("Q", 2.0))
        material = unit.model.Material(name="MAT")
        material.Elastic(table=[[1.0, 0.3]])
        unit.assign_beam_section("STRUTS", "MAT", 0.07, pre_integrate=True)
        self.assertEqual(unit.model.sections["STRUTS"].integration,
                         BEFORE_ANALYSIS)


class TestShellSection(AbaqusTestCase):

    def test_shell_section_is_integrated_during_analysis(self):
        """assign_shell_section integrates through the thickness during the
        analysis (preIntegrate OFF) so nonlinear materials are supported."""
        unit = _build(lambda: fixtures.parallel_plates_3d("SHS", 2.0),
                      element_size=1.0)
        material = unit.model.Material(name="MAT")
        material.Elastic(table=[[1.0, 0.3]])
        unit.assign_shell_section("SHELLS", "MAT", 0.05)
        section = unit.model.sections["SHELLS"]
        self.assertEqual(section.preIntegrate, OFF)
        self.assertAlmostEqual(section.thickness, 0.05)

    def test_shell_section_can_be_pre_integrated(self):
        """pre_integrate=True pre-integrates the section (preIntegrate ON)."""
        unit = _build(lambda: fixtures.parallel_plates_3d("SHP", 2.0),
                      element_size=1.0)
        material = unit.model.Material(name="MAT")
        material.Elastic(table=[[1.0, 0.3]])
        unit.assign_shell_section("SHELLS", "MAT", 0.05, pre_integrate=True)
        self.assertEqual(unit.model.sections["SHELLS"].preIntegrate, ON)


class TestAffineBoundaryConditions(AbaqusTestCase):

    def test_beam_elements_are_rejected(self):
        """Affine displacement BCs do not support rotational DOFs; a beam
        lattice must raise NotImplementedError (periodic BCs handle it)."""
        unit = _build(lambda: fixtures.simple_cubic("AF1", 2.0))
        material = unit.model.Material(name="M")
        material.Elastic(table=[[210.0, 0.3]])
        unit.assign_beam_section("STRUTS", "M", 0.05)
        unit.model.StaticStep(name="LOAD", previous="Initial")
        self.assertRaises(NotImplementedError,
                          unit.apply_affine_displacement_BC,
                          "LOAD", _uniaxial_history())

    def test_warping_elements_are_rejected(self):
        """RVE construction rejects elements with warping DOFs."""
        self.assertRaises(NotImplementedError,
                          _build,
                          lambda: fixtures.simple_cubic("AF2", 2.0),
                          custom_types={1: "B31OS"})

    def test_continuum_is_supported(self):
        """A pure continuum RVE has no rotational DOFs, so the affine BC
        applies the translational Dirichlet BCs without raising."""
        unit = _build(lambda: fixtures.solid_cube("AF3", 1.0))
        material = unit.model.Material(name="M")
        material.Elastic(table=[[100.0, 0.3]])
        unit.assign_solid_section("SOLID", "M")
        unit.model.StaticStep(name="LOAD", previous="Initial")
        unit.apply_affine_displacement_BC("LOAD", _uniaxial_history())
        has_dirichlet = False
        for name in unit.model.boundaryConditions.keys():
            if name.startswith("DIRICHLET_"):
                has_dirichlet = True
                break
        self.assertTrue(has_dirichlet)
