"""
Tests for the automatic boundary detection and the GmshRVE preparation.
These use Gmsh but not the Abaqus kernel.

"""
import unittest

import gmsh

from abqhom.gmshmodel import GmshRVE, Boundary
from abqhom.constants import POINT_GROUPS, EDGE_GROUPS, PERIODIC_PAIRS
from abqhom.tests import fixtures


def _nonempty_groups(boundary, names):
    return [n for n in names if len(boundary[n]) > 0]


class TestBoundaryDetection(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        if not gmsh.is_initialized():
            gmsh.initialize()

    def tearDown(self):
        # remove every model so tests do not interfere
        for name in list(gmsh.model.list()):
            gmsh.model.set_current(name)
            gmsh.model.remove()

    def test_corner_cell_uses_point_groups(self):
        """A cube lattice is periodic and populates the corner groups."""
        fixtures.simple_cubic("TC", length=1.0)
        boundary = Boundary("TC")
        self.assertTrue(boundary.periodic)
        self.assertEqual(len(_nonempty_groups(boundary, POINT_GROUPS)), 8)

    def test_noncorner_cell_uses_edge_groups(self):
        """A cell that crosses the boundary away from the corners is still
        periodic, with empty POINT groups and populated EDGE groups."""
        fixtures.plus_lattice_2d("TP", size=1.0)
        boundary = Boundary("TP")
        self.assertTrue(boundary.periodic)
        self.assertEqual(len(_nonempty_groups(boundary, POINT_GROUPS)), 0)
        self.assertEqual(len(_nonempty_groups(boundary, EDGE_GROUPS[0:4])), 4)

    def test_periodic_pairs_have_equal_counts(self):
        """After meshing, each periodic node pair has matching counts."""
        fixtures.simple_cubic("TM", length=1.0)
        rve = GmshRVE("TM")
        rve.add_periodic_mesh_constraints()
        rve.generate_mesh(element_size=0.5)
        rve.build_node_sets()
        node_sets = rve.node_sets
        for slave, master in PERIODIC_PAIRS.items():
            a, b = node_sets.get(master), node_sets.get(slave)
            if a is not None and b is not None and len(a) and len(b):
                self.assertEqual(len(a), len(b))


class TestGmshRVE(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        if not gmsh.is_initialized():
            gmsh.initialize()

    def tearDown(self):
        for name in list(gmsh.model.list()):
            gmsh.model.set_current(name)
            gmsh.model.remove()

    def test_model_dimension(self):
        """The geometric dimension is detected correctly (2d vs 3d)."""
        fixtures.square_lattice("D2", size=1.0)
        self.assertEqual(GmshRVE("D2").model_dimension, 2)
        fixtures.simple_cubic("D3", length=1.0)
        self.assertEqual(GmshRVE("D3").model_dimension, 3)

    def test_sets_and_mesh(self):
        """Meshing populates the aggregate node and element sets."""
        fixtures.simple_cubic("SM", length=1.0)
        rve = GmshRVE("SM")
        rve.add_periodic_mesh_constraints()
        rve.generate_mesh(element_size=0.5)
        rve.build_element_sets()
        rve.build_node_sets()
        self.assertIn("ALL_NODES", rve.node_sets)
        self.assertIn("ALL_ELEMENTS", rve.element_sets)
        self.assertGreater(len(rve.node_sets["ALL_NODES"]), 0)
        self.assertGreater(len(rve.element_sets["STRUTS"]), 0)

    def test_one_dimensional_model_raises(self):
        """A degenerate (1d) model is rejected by GmshRVE."""
        gmsh.model.add("LINE")
        g = gmsh.model.geo
        line = g.add_line(g.add_point(0.0, 0.0, 0.0),
                          g.add_point(1.0, 0.0, 0.0))
        g.synchronize()
        grp = gmsh.model.add_physical_group(1, [line])
        gmsh.model.set_physical_name(1, grp, "STRUTS")
        self.assertRaises(RuntimeError, GmshRVE, "LINE")


if __name__ == "__main__":
    unittest.main()
