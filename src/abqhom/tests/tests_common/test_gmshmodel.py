"""
Granular tests for individual GmshRVE methods. These cross-check the
abqhom results against Gmsh's own query functions (mesh nodes/elements,
periodic-node correspondences, entity coordinates) rather than testing a
whole workflow. Gmsh only, no Abaqus kernel.

"""
import os
import unittest

import numpy as np

import gmsh

from abqhom.constants import (POINT_GROUPS, EDGE_GROUPS, FACE_GROUPS,
                              PERIODIC_PAIRS)
from abqhom.gmshmodel import GmshRVE
from abqhom.tests import fixtures


class _GmshTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        if not gmsh.is_initialized():
            gmsh.initialize()

    def tearDown(self):
        for name in list(gmsh.model.list()):
            gmsh.model.set_current(name)
            gmsh.model.remove()


class TestDimensions(_GmshTestCase):

    def test_model_dimension(self):
        fixtures.square_lattice("MD2", size=1.0)
        self.assertEqual(GmshRVE("MD2").model_dimension, 2)
        fixtures.simple_cubic("MD3", length=1.0)
        self.assertEqual(GmshRVE("MD3").model_dimension, 3)

    def test_mesh_dimension_is_one_for_struts(self):
        fixtures.simple_cubic("MESHDIM", length=2.0)
        rve = GmshRVE("MESHDIM")
        rve.generate_mesh(element_size=1.0)
        self.assertEqual(rve.mesh_dimension, 1)


class TestNodesAndElements(_GmshTestCase):

    def test_get_nodes_are_unique_and_include_corners(self):
        fixtures.simple_cubic("NODES", length=2.0)
        rve = GmshRVE("NODES")
        rve.generate_mesh(element_size=1.0)
        tags, coords = rve.get_nodes()

        # No duplicate node tags and matching coordinate array shape
        self.assertEqual(len(tags), len(set(tags.tolist())))
        self.assertEqual(coords.shape, (len(tags), 3))

        # The eight cube corners must be among the returned coordinates
        corners = [(x, y, z) for x in (-1.0, 1.0)
                   for y in (-1.0, 1.0) for z in (-1.0, 1.0)]
        for corner in corners:
            distances = np.linalg.norm(coords - np.array(corner), axis=1)
            self.assertTrue(np.min(distances) < 1e-9)

    def test_get_elements_are_two_node_lines(self):
        fixtures.simple_cubic("ELEMS", length=2.0)
        rve = GmshRVE("ELEMS")
        rve.generate_mesh(element_size=2.0)  # one element per strut
        element_tags, connectivity = rve.get_elements()

        # Gmsh element type 1 is the 2-node line
        self.assertIn(1, element_tags)
        self.assertEqual(connectivity[1].shape[1], 2)
        self.assertEqual(len(element_tags[1]), 12)  # twelve cube edges

    def test_second_order_mesh_gives_three_node_lines(self):
        fixtures.simple_cubic("ORDER2", length=2.0)
        rve = GmshRVE("ORDER2")
        rve.generate_mesh(element_order=2, element_size=2.0)
        element_tags, connectivity = rve.get_elements()

        # Gmsh element type 8 is the 3-node second-order line
        self.assertIn(8, element_tags)
        self.assertEqual(connectivity[8].shape[1], 3)


class TestPeriodicity(_GmshTestCase):

    def test_periodic_node_sets_are_offset_by_the_lattice_vector(self):
        """A master node and its paired slave node differ by the lattice
        vector, verified against the actual Gmsh node coordinates."""
        length = 2.0
        fixtures.simple_cubic("PNS", length=length)
        rve = GmshRVE("PNS")
        rve.add_periodic_mesh_constraints()
        rve.generate_mesh(element_size=1.0)
        rve.build_node_sets()

        master = rve.node_sets["POINT_1"]   # (xmin, ymin, zmin)
        slave = rve.node_sets["POINT_7"]    # (xmax, ymax, zmax)
        self.assertEqual(len(master), 1)
        self.assertEqual(len(slave), 1)

        master_coord = gmsh.model.mesh.get_node(int(master[0]))[0]
        slave_coord = gmsh.model.mesh.get_node(int(slave[0]))[0]
        offset = np.array(slave_coord) - np.array(master_coord)
        self.assertTrue(np.allclose(offset, [length, length, length]))

    def test_periodic_mesh_constraint_is_registered_in_gmsh(self):
        """add_periodic_mesh_constraints actually sets a Gmsh periodic
        correspondence on the slave faces of a continuum cube."""
        fixtures.solid_cube("PMC", length=2.0)
        rve = GmshRVE("PMC")
        rve.add_periodic_mesh_constraints()
        rve.generate_mesh(element_size=1.0)
        self.assertTrue(rve.periodic_mesh)

        # FACE_2 (z = zmax) is the slave of FACE_1 (z = zmin)
        slave_entities = rve.boundary["FACE_2"]
        self.assertTrue(len(slave_entities) > 0)
        dim, tag = slave_entities[0]
        _, node_tags, _, _ = gmsh.model.mesh.get_periodic_nodes(dim, tag)
        self.assertGreater(len(node_tags), 0)


class TestNodeSets(_GmshTestCase):

    def test_boundary_and_inner_nodes_partition_all_nodes(self):
        fixtures.simple_cubic("PART", length=2.0)
        rve = GmshRVE("PART")
        rve.add_periodic_mesh_constraints()
        rve.generate_mesh(element_size=1.0)
        rve.build_node_sets()

        all_nodes = set(rve.node_sets["ALL_NODES"].tolist())
        boundary = set(rve.node_sets["BOUNDARY_NODES"].tolist())
        inner = set(rve.node_sets["INNER_NODES"].tolist())
        self.assertEqual(boundary | inner, all_nodes)
        self.assertEqual(boundary & inner, set())

    def test_element_sets_follow_physical_groups(self):
        fixtures.bcc("ESETS", length=2.0)
        rve = GmshRVE("ESETS")
        rve.generate_mesh(element_size=3.0)  # one element per strut
        rve.build_element_sets()
        self.assertIn("CUBIC", rve.element_sets)
        self.assertIn("BCC", rve.element_sets)
        self.assertEqual(len(rve.element_sets["CUBIC"]), 12)  # twelve edges
        self.assertEqual(len(rve.element_sets["BCC"]), 8)     # eight diagonals

    def test_boundary_element_sets_classify_cubic_edges(self):
        fixtures.bcc("BOUNDARY_ESETS", length=2.0)
        rve = GmshRVE("BOUNDARY_ESETS")
        rve.generate_mesh(element_size=1.0)  # two elements per strut
        rve.build_element_sets()

        edge_elements = set()
        for group in EDGE_GROUPS:
            self.assertEqual(len(rve.element_sets[group]), 2)
            edge_elements.update(rve.element_sets[group].tolist())

        self.assertEqual(edge_elements,
                         set(rve.element_sets["CUBIC"].tolist()))
        self.assertTrue(edge_elements.isdisjoint(
            set(rve.element_sets["BCC"].tolist())))
        for group in POINT_GROUPS + FACE_GROUPS:
            self.assertEqual(len(rve.element_sets[group]), 0)

        self.assertEqual(set(rve.element_sets["EDGE_ELEMENTS"]),
                         edge_elements)
        self.assertEqual(set(rve.element_sets["BOUNDARY_ELEMENTS"]),
                         edge_elements)
        self.assertEqual(set(rve.element_sets["INNER_ELEMENTS"]),
                         set(rve.element_sets["BCC"]))
        self.assertEqual(len(rve.element_sets["POINT_ELEMENTS"]), 0)
        self.assertEqual(len(rve.element_sets["FACE_ELEMENTS"]), 0)

        secondary_elements = set()
        for group in PERIODIC_PAIRS:
            secondary_elements.update(rve.element_sets[group])
        self.assertEqual(set(rve.element_sets["SECONDARY_ELEMENTS"]),
                         secondary_elements)

    def test_boundary_element_sets_support_two_dimensions(self):
        fixtures.square_lattice("BOUNDARY_ESETS_2D", size=2.0)
        rve = GmshRVE("BOUNDARY_ESETS_2D")
        rve.generate_mesh(element_size=1.0)
        rve.build_element_sets()

        boundary_elements = set()
        for group in EDGE_GROUPS[:4]:
            self.assertEqual(len(rve.element_sets[group]), 2)
            boundary_elements.update(rve.element_sets[group].tolist())

        self.assertEqual(boundary_elements,
                         set(rve.element_sets["STRUTS"].tolist()))
        for group in POINT_GROUPS + EDGE_GROUPS[4:] + FACE_GROUPS:
            self.assertEqual(len(rve.element_sets[group]), 0)

        secondary_elements = (set(rve.element_sets["EDGE_2"])
                              | set(rve.element_sets["EDGE_3"]))
        self.assertEqual(set(rve.element_sets["SECONDARY_ELEMENTS"]),
                         secondary_elements)

    def test_continuum_elements_are_not_boundary_elements(self):
        fixtures.solid_cube("CONTINUUM_ESETS", length=2.0)
        rve = GmshRVE("CONTINUUM_ESETS")
        rve.generate_mesh(element_size=1.0)
        rve.build_element_sets()

        for group in POINT_GROUPS + EDGE_GROUPS + FACE_GROUPS:
            self.assertEqual(len(rve.element_sets[group]), 0)
        self.assertEqual(len(rve.element_sets["BOUNDARY_ELEMENTS"]), 0)
        self.assertEqual(set(rve.element_sets["INNER_ELEMENTS"]),
                         set(rve.element_sets["ALL_ELEMENTS"]))


class TestWriteAbaqusInput(_GmshTestCase):

    def test_input_file_contains_nodes_and_beam_elements(self):
        fixtures.simple_cubic("WINP", length=2.0)
        rve = GmshRVE("WINP")
        rve.generate_mesh(element_size=2.0)
        path = rve.write_abaqus_input()
        try:
            content = open(path, "r").read()
        finally:
            if os.path.isfile(path):
                os.remove(path)
        self.assertIn("*NODE", content)
        # Gmsh 2-node line maps to the Abaqus B33 beam in 3d
        self.assertIn("*ELEMENT, TYPE=B33", content)

    def test_custom_types_do_not_mutate_the_shared_type_map(self):
        """Passing custom_types must not mutate the shared default element
        type map (regression guard for the in-place update bug)."""
        from abqhom.constants import ABAQUS_ELEMENT_TYPES_3D
        before = dict(ABAQUS_ELEMENT_TYPES_3D)
        fixtures.simple_cubic("WINP2", length=2.0)
        rve = GmshRVE("WINP2")
        rve.generate_mesh(element_size=2.0)
        path = rve.write_abaqus_input(custom_types={1: "B31"})
        try:
            content = open(path, "r").read()
        finally:
            if os.path.isfile(path):
                os.remove(path)
        # The override reaches the written file ...
        self.assertIn("*ELEMENT, TYPE=B31", content)
        # ... but the shared default map is left untouched.
        self.assertEqual(ABAQUS_ELEMENT_TYPES_3D, before)

    def test_unsupported_element_type_is_rejected(self):
        """Mapping an element to an unsupported type (a 5-DOF thin shell)
        raises rather than writing an unusable input file."""
        fixtures.simple_cubic("WINP3", length=2.0)
        rve = GmshRVE("WINP3")
        rve.generate_mesh(element_size=2.0)
        self.assertRaises(NotImplementedError, rve.write_abaqus_input,
                          {1: "S8R5"})


class TestReadFromMesh(_GmshTestCase):

    def test_pairs_periodic_faces_of_a_prebuilt_mesh(self):
        """read_from_mesh classifies the nodes of a ready-made periodic mesh
        onto the box groups, detects periodicity and pairs opposite faces
        one to one."""
        L = 2.0
        fixtures.discrete_periodic_plates("DPP", size=L)
        rve = GmshRVE.read_from_mesh("DPP")

        self.assertEqual(rve.model_dimension, 3)
        self.assertTrue(rve.periodic_mesh)
        ns = rve.node_sets
        self.assertEqual(len(ns["ALL_NODES"]), 18)
        self.assertIn("PLATES", ns)

        # The two plates lie on the z-faces: FACE_1 (z=0) and FACE_2 (z=L)
        # each hold the single interior node and must pair up.
        self.assertEqual(len(ns["FACE_1"]), len(ns["FACE_2"]))
        self.assertGreater(len(ns["FACE_1"]), 0)

        # The shell elements are classified on their respective faces.
        es = rve.element_sets
        self.assertEqual(len(es["FACE_1"]), 8)
        self.assertEqual(len(es["FACE_2"]), 8)
        face_elements = set(es["FACE_1"]) | set(es["FACE_2"])
        self.assertEqual(face_elements, set(es["PLATES"]))

        # The paired nodes differ by exactly the lattice vector in z.
        tags, coords, _ = gmsh.model.mesh.get_nodes()
        coords = np.array(coords).reshape(-1, 3)
        coord = {int(t): coords[i] for i, t in enumerate(tags)}
        for master_tag, slave_tag in zip(ns["FACE_1"], ns["FACE_2"]):
            offset = coord[int(slave_tag)] - coord[int(master_tag)]
            self.assertTrue(np.allclose(offset, [0.0, 0.0, L]))

    def test_detects_a_non_periodic_mesh(self):
        """A mesh whose opposite faces are not node-matched is reported as
        non-periodic (periodic_mesh is False)."""
        fixtures.discrete_periodic_plates("DPP2", size=2.0, perturb=0.3)
        rve = GmshRVE.read_from_mesh("DPP2")
        self.assertEqual(rve.model_dimension, 3)
        self.assertFalse(rve.periodic_mesh)


if __name__ == "__main__":
    unittest.main()
