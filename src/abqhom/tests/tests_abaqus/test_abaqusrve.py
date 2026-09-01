"""
Granular tests for AbaqusRVE that exercise individual pieces of the API
without running a finite element analysis: construction, set transfer,
the periodic sharing weights, the rigid-body anchor selection and the
recorded section radii. These need the Abaqus kernel.

"""
import unittest

import numpy as np

import gmsh

from abqhom.utils import voigt_to_tensor
from abqhom.gmshmodel import GmshRVE
from abqhom.abqmodel import AbaqusRVE, LoadHistory
from abqhom.tests import fixtures

from abaqusConstants import (THREE_D, TWO_D_PLANAR, ON, OFF, DURING_ANALYSIS,
                             BEFORE_ANALYSIS)


def _uniaxial_history():
    """A trivial uniaxial (symmetric) strain history for BC tests."""
    def strain_history(time):
        strain = np.zeros(6)
        strain[0] = 1.0e-3
        return time*voigt_to_tensor(strain, "STRAIN")
    return LoadHistory(np.array([0.0, 1.0]), strain_history)


def _build(builder, element_size=1.0, custom_types=None):
    name = builder()
    rve = GmshRVE(name)
    rve.add_periodic_mesh_constraints()
    rve.generate_mesh(element_size=element_size)
    return AbaqusRVE(rve, custom_types or {})


class _AbaqusTestCase(unittest.TestCase):

    def tearDown(self):
        for name in list(gmsh.model.list()):
            gmsh.model.set_current(name)
            gmsh.model.remove()


class TestConstruction(_AbaqusTestCase):

    def test_model_dimension_3d(self):
        unit = _build(lambda: fixtures.simple_cubic("A", 2.0))
        self.assertEqual(unit.model_dimension, THREE_D)
        self.assertIsNotNone(unit.part)
        self.assertIsNotNone(unit.instance)

    def test_model_dimension_2d(self):
        unit = _build(lambda: fixtures.square_lattice("B", 2.0))
        self.assertEqual(unit.model_dimension, TWO_D_PLANAR)

    def test_bounding_box_matches_geometry(self):
        unit = _build(lambda: fixtures.simple_cubic("C", 2.0))
        self.assertTrue(np.allclose(unit.get_bounding_box(),
                                    [-1, -1, -1, 1, 1, 1]))

    def test_node_and_element_sets_are_transferred(self):
        unit = _build(lambda: fixtures.simple_cubic("D", 2.0))
        self.assertIn("ELSET_STRUTS", unit.part.sets)
        self.assertIn("NSET_ALL_NODES", unit.part.sets)
        self.assertIn("NSET_POINT_1", unit.part.sets)


class TestSharingWeights(_AbaqusTestCase):

    def test_cube_edge_struts_get_one_quarter(self):
        """Every strut of the simple-cubic cell lies on a bounding-box
        edge (two faces) and is therefore weighted 1/4."""
        unit = _build(lambda: fixtures.simple_cubic("E", 2.0))
        weights = unit._element_sharing_weights()
        self.assertEqual(set(weights.values()), set([0.25]))

    def test_square_edge_struts_get_one_half(self):
        """Each strut of the 2d square cell lies on one bounding-box face
        and is weighted 1/2."""
        unit = _build(lambda: fixtures.square_lattice("F", 2.0))
        weights = unit._element_sharing_weights()
        self.assertEqual(set(weights.values()), set([0.5]))

    def test_transversal_struts_are_not_weighted(self):
        """Struts that only cross the boundary keep weight 1."""
        unit = _build(lambda: fixtures.plus_lattice_2d("G", 2.0))
        weights = unit._element_sharing_weights()
        self.assertEqual(set(weights.values()), set([1.0]))

    def test_bcc_has_mixed_weights(self):
        """The BCC cell mixes edge struts (1/4) and interior diagonals
        (1)."""
        unit = _build(lambda: fixtures.bcc("H", 2.0))
        values = list(unit._element_sharing_weights().values())
        self.assertIn(0.25, values)
        self.assertIn(1.0, values)
        self.assertNotIn(0.5, values)


class TestRigidBodyAnchor(_AbaqusTestCase):

    def test_anchor_is_a_corner_for_a_cubic_cell(self):
        unit = _build(lambda: fixtures.simple_cubic("K", 2.0))
        region, node = unit._select_rigid_body_anchor()
        corner_labels = [n.label for n
                         in unit.instance.sets["NSET_POINT_1"].nodes]
        self.assertIn(node.label, corner_labels)

    def test_anchor_is_not_a_periodic_slave(self):
        """The auto-selected anchor must be a pure-master node so it can
        carry a Dirichlet BC without over-constraining a periodic
        equation."""
        unit = _build(lambda: fixtures.plus_lattice_2d("L", 2.0))
        _, node = unit._select_rigid_body_anchor()
        from abqhom.constants import PERIODIC_PAIRS
        slave_labels = set()
        for slave in PERIODIC_PAIRS.keys():
            set_name = "NSET_" + slave
            if set_name in unit.instance.sets:
                for n in unit.instance.sets[set_name].nodes:
                    slave_labels.add(n.label)
        self.assertNotIn(node.label, slave_labels)

    def test_user_defined_anchor_is_used(self):
        # The BCC cell has an interior center node, so INNER_NODES is not
        # empty and can be used as a user-defined anchor.
        unit = _build(lambda: fixtures.bcc("M", 2.0))
        _, node = unit._select_rigid_body_anchor(fixed_node="INNER_NODES")
        inner_labels = [n.label for n
                        in unit.instance.sets["NSET_INNER_NODES"].nodes]
        self.assertIn(node.label, inner_labels)


class TestSectionRadius(_AbaqusTestCase):

    def test_beam_section_records_the_radius(self):
        unit = _build(lambda: fixtures.simple_cubic("N", 2.0))
        material = unit.model.Material(name="MAT")
        material.Elastic(table=[[1.0, 0.3]])
        unit.assign_beam_section("STRUTS", "MAT", 0.07)
        self.assertAlmostEqual(unit._section_radius["STRUTS"], 0.07)

    def test_truss_section_records_the_equivalent_radius(self):
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


class TestShellSection(_AbaqusTestCase):

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


class TestAffineBoundaryConditions(_AbaqusTestCase):

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
        """Open-section beams (rotational and warping DOFs) are likewise
        rejected by the affine BC."""
        unit = _build(lambda: fixtures.simple_cubic("AF2", 2.0),
                      custom_types={1: "B31OS"})
        unit.model.StaticStep(name="LOAD", previous="Initial")
        self.assertRaises(NotImplementedError,
                          unit.apply_affine_displacement_BC,
                          "LOAD", _uniaxial_history())

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


if __name__ == "__main__":
    unittest.main()
