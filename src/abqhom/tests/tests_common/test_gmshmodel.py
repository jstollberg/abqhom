# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Granular tests for individual GmshRVE methods. These cross-check the
abqhom results against Gmsh's own query functions (mesh nodes/elements,
periodic-node correspondences, entity coordinates) rather than testing a
whole workflow. Gmsh only, no Abaqus kernel.

"""
import os

import numpy as np

import gmsh

from abqhom.constants import (POINT_GROUPS, EDGE_GROUPS, FACE_GROUPS,
                              PERIODIC_PAIRS, ABAQUS_CONNECTIVITY_ORDER)
from abqhom.gmshmodel import Boundary, GmshRVE
from abqhom.tests import fixtures
from abqhom.tests.base import GmshTestCase


class TestDimensions(GmshTestCase):
    """Test geometric and mesh dimensions of Gmsh RVEs."""

    def test_boundary_must_belong_to_rve_model(self):
        """A supplied boundary must describe the requested Gmsh model."""
        fixtures.simple_cubic("BOUNDARY_MODEL", length=1.0)
        boundary = Boundary("BOUNDARY_MODEL")
        fixtures.simple_cubic("RVE_MODEL", length=1.0)
        self.assertRaises(ValueError, GmshRVE, "RVE_MODEL", boundary)

    def test_model_dimension(self):
        """Planar and spatial models report dimensions two and three."""
        fixtures.square_lattice("MD2", size=1.0)
        dimension_2d = GmshRVE("MD2").model_dimension
        self.assertEqual(dimension_2d, 2)
        self.assertIs(type(dimension_2d), int)
        fixtures.simple_cubic("MD3", length=1.0)
        dimension_3d = GmshRVE("MD3").model_dimension
        self.assertEqual(dimension_3d, 3)
        self.assertIs(type(dimension_3d), int)

    def test_mesh_dimension_is_one_for_struts(self):
        """A strut-only cell reports a one-dimensional mesh."""
        fixtures.simple_cubic("MESHDIM", length=2.0)
        rve = GmshRVE("MESHDIM")
        rve.generate_mesh(element_size=1.0)
        self.assertEqual(rve.mesh_dimension, 1)

    def test_string_representation_contains_model_name(self):
        """The string representation identifies the underlying model."""
        fixtures.square_lattice("REPR", size=1.0)
        self.assertIn("GmshRVE 'REPR'", str(GmshRVE("REPR")))


class TestNodesAndElements(GmshTestCase):
    """Test retrieval of nodes, elements, and connectivity."""

    def test_get_nodes_are_unique_and_include_corners(self):
        """Node extraction is unique and includes every cell corner."""
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
        """A first-order strut mesh contains two-node line elements."""
        fixtures.simple_cubic("ELEMS", length=2.0)
        rve = GmshRVE("ELEMS")
        rve.generate_mesh(element_size=2.0)  # one element per strut
        element_tags, connectivity = rve.get_elements()

        # Gmsh element type 1 is the 2-node line
        self.assertIn(1, element_tags)
        self.assertEqual(connectivity[1].shape[1], 2)
        self.assertEqual(len(element_tags[1]), 12)  # twelve cube edges

    def test_second_order_mesh_gives_three_node_lines(self):
        """A second-order strut mesh contains three-node line elements."""
        fixtures.simple_cubic("ORDER2", length=2.0)
        rve = GmshRVE("ORDER2")
        rve.generate_mesh(element_order=2, element_size=2.0)
        element_tags, connectivity = rve.get_elements()

        # Gmsh element type 8 is the 3-node second-order line
        self.assertIn(8, element_tags)
        self.assertEqual(connectivity[8].shape[1], 3)

    def test_missing_connectivity_mapping_raises_clear_error(self):
        """Unsupported Gmsh connectivity raises a descriptive exception."""
        fixtures.simple_cubic("NO_CONNECTIVITY", length=2.0)
        rve = GmshRVE("NO_CONNECTIVITY")
        rve.generate_mesh(element_size=2.0)
        connectivity_order = ABAQUS_CONNECTIVITY_ORDER.pop(1)
        try:
            self.assertRaises(NotImplementedError, rve.get_elements)
        finally:
            ABAQUS_CONNECTIVITY_ORDER[1] = connectivity_order

    def test_mesh_operations_activate_their_rve_model(self):
        """Meshing and element access do not act on another current model."""
        fixtures.simple_cubic("ACTIVE_RVE", length=2.0)
        rve = GmshRVE("ACTIVE_RVE")
        fixtures.simple_cubic("OTHER_MODEL", length=2.0)

        rve.generate_mesh(element_size=2.0)
        self.assertEqual(gmsh.model.get_current(), "ACTIVE_RVE")

        gmsh.model.set_current("OTHER_MODEL")
        element_tags, _ = rve.get_elements()
        self.assertEqual(gmsh.model.get_current(), "ACTIVE_RVE")
        self.assertIn(1, element_tags)

    def test_generate_mesh_restores_global_options(self):
        """Per-call meshing settings do not leak into later Gmsh models."""
        fixtures.simple_cubic("RESTORE_OPTIONS", length=2.0)
        rve = GmshRVE("RESTORE_OPTIONS")
        option_values = {
            "Mesh.ElementOrder": 1.0,
            "Mesh.SecondOrderIncomplete": 0.0,
            "Mesh.Algorithm": 6.0,
            "Mesh.Algorithm3D": 1.0,
            "Mesh.MeshSizeFromPoints": 1.0,
            "Mesh.MeshSizeMin": 0.125,
            "Mesh.MeshSizeMax": 4.0,
        }
        for name, value in option_values.items():
            gmsh.option.set_number(name, value)

        rve.generate_mesh(element_order=2, element_size=0.75,
                          algorithm_2d=5, algorithm_3d=4)

        for name, expected in option_values.items():
            self.assertEqual(gmsh.option.get_number(name), expected)


class TestPeriodicity(GmshTestCase):
    """Test periodic mesh constraints and node correspondence."""

    def test_periodic_node_sets_are_offset_by_the_lattice_vector(self):
        """Paired primary and secondary nodes differ by a lattice vector."""
        length = 2.0
        fixtures.simple_cubic("PNS", length=length)
        rve = GmshRVE("PNS")
        rve.add_periodic_mesh_constraints()
        rve.generate_mesh(element_size=1.0)
        rve.build_node_sets()

        primary = rve.node_sets["POINT_1"]    # (xmin, ymin, zmin)
        secondary = rve.node_sets["POINT_7"]  # (xmax, ymax, zmax)
        self.assertEqual(len(primary), 1)
        self.assertEqual(len(secondary), 1)

        primary_coord = gmsh.model.mesh.get_node(int(primary[0]))[0]
        secondary_coord = gmsh.model.mesh.get_node(int(secondary[0]))[0]
        offset = np.array(secondary_coord) - np.array(primary_coord)
        self.assertTrue(np.allclose(offset, [length, length, length]))

    def test_periodic_mesh_constraint_is_registered_in_gmsh(self):
        """Gmsh records periodic correspondence on secondary faces."""
        fixtures.solid_cube("PMC", length=2.0)
        rve = GmshRVE("PMC")
        rve.add_periodic_mesh_constraints()
        rve.generate_mesh(element_size=1.0)
        self.assertTrue(rve.periodic_mesh)

        # FACE_2 (z = zmax) is secondary to FACE_1 (z = zmin)
        secondary_entities = rve.boundary["FACE_2"]
        self.assertTrue(len(secondary_entities) > 0)
        dim, tag = secondary_entities[0]
        _, node_tags, _, _ = gmsh.model.mesh.get_periodic_nodes(dim, tag)
        self.assertGreater(len(node_tags), 0)

    def test_periodic_node_matching_preserves_multiplicity(self):
        """Periodic groups with unequal node counts must not match."""
        signatures = {"POINT_1": (-1, -1, -1),
                      "POINT_2": (1, -1, -1)}
        coordinates = {"POINT_1": np.empty((1, 0)),
                       "POINT_2": np.empty((2, 0))}
        self.assertFalse(GmshRVE._periodic_node_sets_match(
            coordinates, signatures, tolerance=1.0e-12))

    def test_periodic_node_matching_detects_small_scale_mismatch(self):
        """Matching does not erase differences through decimal rounding."""
        signatures = {"EDGE_1": (0, -1, -1),
                      "EDGE_3": (0, 1, -1)}
        coordinates = {
            "EDGE_1": np.array([[0.0], [2.0e-10]]),
            "EDGE_3": np.array([[0.0], [3.0e-10]]),
        }
        self.assertFalse(GmshRVE._periodic_node_sets_match(
            coordinates, signatures, tolerance=1.0e-15))


class TestNodeAndElementSets(GmshTestCase):
    """Test generated node and element sets."""

    def test_boundary_and_inner_nodes_partition_all_nodes(self):
        """Boundary and inner node sets form a disjoint partition."""
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
        """Element sets preserve the Gmsh physical-group membership."""
        fixtures.bcc("ESETS", length=2.0)
        rve = GmshRVE("ESETS")
        rve.generate_mesh(element_size=3.0)  # one element per strut
        rve.build_element_sets()
        self.assertIn("CUBIC", rve.element_sets)
        self.assertIn("BCC", rve.element_sets)
        self.assertEqual(len(rve.element_sets["CUBIC"]), 12)  # twelve edges
        self.assertEqual(len(rve.element_sets["BCC"]), 8)     # eight diagonals

    def test_boundary_element_sets_classify_cubic_edges(self):
        """Spatial boundary struts are classified by their cell edges."""
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
        """Planar boundary struts are classified by their cell edges."""
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
        """Volume elements remain inner despite touching the boundary."""
        fixtures.solid_cube("CONTINUUM_ESETS", length=2.0)
        rve = GmshRVE("CONTINUUM_ESETS")
        rve.generate_mesh(element_size=1.0)
        rve.build_element_sets()

        for group in POINT_GROUPS + EDGE_GROUPS + FACE_GROUPS:
            self.assertEqual(len(rve.element_sets[group]), 0)
        self.assertEqual(len(rve.element_sets["BOUNDARY_ELEMENTS"]), 0)
        self.assertEqual(set(rve.element_sets["INNER_ELEMENTS"]),
                         set(rve.element_sets["ALL_ELEMENTS"]))


class TestWriteAbaqusInput(GmshTestCase):
    """Test Abaqus input-file generation and type validation."""

    def test_input_file_contains_nodes_and_beam_elements(self):
        """Abaqus export writes nodes and the default beam element type."""
        fixtures.simple_cubic("WINP", length=2.0)
        rve = GmshRVE("WINP")
        rve.generate_mesh(element_size=2.0)
        path = rve.write_abaqus_input()
        try:
            with open(path, "r") as input_file:
                content = input_file.read()
        finally:
            if os.path.isfile(path):
                os.remove(path)
        self.assertIn("*NODE", content)

        # Gmsh 2-node line maps to the Abaqus B33 beam in 3D
        self.assertIn("*ELEMENT, TYPE=B33", content)

    def test_custom_types_do_not_mutate_the_shared_type_map(self):
        """Custom types do not mutate the shared default type map."""
        from abqhom.constants import ABAQUS_ELEMENT_TYPES_3D
        before = dict(ABAQUS_ELEMENT_TYPES_3D)
        fixtures.simple_cubic("WINP2", length=2.0)
        rve = GmshRVE("WINP2")
        rve.generate_mesh(element_size=2.0)
        path = rve.write_abaqus_input(custom_types={1: "B31"})
        try:
            with open(path, "r") as input_file:
                content = input_file.read()
        finally:
            if os.path.isfile(path):
                os.remove(path)

        # The override reaches the written file ...
        self.assertIn("*ELEMENT, TYPE=B31", content)

        # ... but the shared default map is left untouched.
        self.assertEqual(ABAQUS_ELEMENT_TYPES_3D, before)

    def test_unsupported_element_type_is_rejected(self):
        """A mapped unsupported element type is rejected."""
        fixtures.simple_cubic("WINP3", length=2.0)
        rve = GmshRVE("WINP3")
        rve.generate_mesh(element_size=2.0)
        self.assertRaises(NotImplementedError, rve.write_abaqus_input,
                          {1: "S8R5"})

    def test_warping_element_type_is_rejected(self):
        """Element types carrying a warping DOF are not supported."""
        fixtures.simple_cubic("WINP4", length=2.0)
        rve = GmshRVE("WINP4")
        rve.generate_mesh(element_size=2.0)
        self.assertRaises(NotImplementedError, rve.write_abaqus_input,
                          {1: "B31OS"})


class TestReadFromMesh(GmshTestCase):
    """Test construction from existing periodic meshes."""

    def test_pairs_periodic_edges_of_a_prebuilt_2d_mesh(self):
        """A periodic planar mesh has matching opposite-edge nodes."""
        size = 2.0
        fixtures.discrete_periodic_square("DPS", size=size)
        rve = GmshRVE.read_from_mesh("DPS")

        self.assertEqual(rve.model_dimension, 2)
        self.assertTrue(rve.periodic_mesh)
        ns = rve.node_sets
        self.assertEqual(len(ns["ALL_NODES"]), 9)
        self.assertIn("SOLID", ns)

        self.assertEqual(len(ns["EDGE_1"]), 1)
        self.assertEqual(len(ns["EDGE_3"]), 1)
        self.assertEqual(len(ns["EDGE_2"]), 1)
        self.assertEqual(len(ns["EDGE_4"]), 1)

        tags, coords, _ = gmsh.model.mesh.get_nodes()
        coords = np.array(coords).reshape(-1, 3)
        coord = {int(tag): coords[i] for i, tag in enumerate(tags)}
        bottom = coord[int(ns["EDGE_1"][0])]
        top = coord[int(ns["EDGE_3"][0])]
        left = coord[int(ns["EDGE_4"][0])]
        right = coord[int(ns["EDGE_2"][0])]
        self.assertTrue(np.allclose(top - bottom, [0.0, size, 0.0]))
        self.assertTrue(np.allclose(right - left, [size, 0.0, 0.0]))

        es = rve.element_sets
        self.assertEqual(set(es["SOLID"]), set(es["ALL_ELEMENTS"]))
        self.assertEqual(len(es["BOUNDARY_ELEMENTS"]), 0)

    def test_detects_a_non_periodic_prebuilt_2d_mesh(self):
        """Mismatched planar boundary nodes are detected as nonperiodic."""
        fixtures.discrete_periodic_square("DPS2", size=2.0, perturb=0.2)
        rve = GmshRVE.read_from_mesh("DPS2")
        self.assertEqual(rve.model_dimension, 2)
        self.assertFalse(rve.periodic_mesh)

    def test_rejects_a_prebuilt_2d_mesh_outside_xy_plane(self):
        """A pre-built x-z mesh is rejected instead of being projected."""
        fixtures.square_lattice_in_plane(
            "DPS_XZ", size=2.0, axes=(0, 2))
        gmsh.model.mesh.generate(1)

        with self.assertRaises(RuntimeError) as context:
            GmshRVE.read_from_mesh("DPS_XZ")
        self.assertIn("x-y plane", str(context.exception))

    def test_pairs_periodic_faces_of_a_prebuilt_mesh(self):
        """A periodic spatial mesh has matching opposite-face nodes."""
        size = 2.0
        fixtures.discrete_periodic_plates("DPP", size=size)
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
        for primary_tag, secondary_tag in zip(ns["FACE_1"], ns["FACE_2"]):
            offset = coord[int(secondary_tag)] - coord[int(primary_tag)]
            self.assertTrue(np.allclose(offset, [0.0, 0.0, size]))

    def test_detects_a_non_periodic_mesh(self):
        """Mismatched spatial boundary nodes are detected as nonperiodic."""
        fixtures.discrete_periodic_plates("DPP2", size=2.0, perturb=0.3)
        rve = GmshRVE.read_from_mesh("DPP2")
        self.assertEqual(rve.model_dimension, 3)
        self.assertFalse(rve.periodic_mesh)
