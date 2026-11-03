# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Granular tests for the Boundary class and the GmshRVE input validation.
Gmsh only, no Abaqus kernel.

"""
import numpy as np

import gmsh

from abqhom.gmshmodel import GmshRVE, Boundary
from abqhom.constants import (POINT_GROUPS, EDGE_GROUPS, FACE_GROUPS,
                              PERIODIC_PAIRS)
from abqhom.tests import fixtures
from abqhom.tests.base import GmshTestCase


def _nonempty(boundary, names):
    """Return the supplied boundary groups that contain entities."""
    return [name for name in names if len(boundary[name]) > 0]


class TestBoundaryGeometry(GmshTestCase):

    def test_bounding_box_from_points(self):
        """Boundary points recover the exact cuboidal bounding box."""
        fixtures.simple_cubic("BB", length=2.0)
        box = Boundary("BB").get_bounding_box_from_points()
        self.assertTrue(np.allclose(box, [-1, -1, -1, 1, 1, 1]))

    def test_isolated_points_do_not_enlarge_cell_bounds(self):
        """Disconnected construction points do not define RVE bounds."""
        fixtures.simple_cubic("ISOLATED_POINT", length=2.0)
        gmsh.model.geo.add_point(10.0, 0.0, 0.0)
        gmsh.model.geo.synchronize()

        boundary = Boundary("ISOLATED_POINT")

        self.assertTrue(boundary.periodic)
        self.assertTrue(np.allclose(
            boundary.get_bounding_box_from_points(),
            [-1, -1, -1, 1, 1, 1]))

    def test_group_coordinate_identifiers(self):
        """Group identifiers encode fixed and unconstrained coordinates."""
        fixtures.simple_cubic("ID", length=2.0)
        identifiers = Boundary("ID").get_group_coordinate_identifiers()

        self.assertEqual(list(identifiers.keys()),
                         POINT_GROUPS + EDGE_GROUPS + FACE_GROUPS)

        # POINT_1 is the minimum corner and POINT_7 the maximum corner.
        self.assertTrue(np.allclose(identifiers["POINT_1"], [-1, -1, -1]))
        self.assertTrue(np.allclose(identifiers["POINT_7"], [1, 1, 1]))

        # FACE_1 lies at z_min and is free in x and y.
        face1 = identifiers["FACE_1"]
        self.assertTrue(np.isnan(face1[0]) and np.isnan(face1[1]))
        self.assertAlmostEqual(face1[2], -1.0)

    def test_static_group_identifiers_match_model_bounds(self):
        """Static and model-based group identifiers are equivalent."""
        fixtures.simple_cubic("STATIC_ID", length=2.0)
        model_identifiers = Boundary(
            "STATIC_ID").get_group_coordinate_identifiers()
        static_identifiers = Boundary.group_coordinate_identifiers(
            [-1.0, -1.0, -1.0], [1.0, 1.0, 1.0])

        for name in model_identifiers:
            self.assertTrue(np.array_equal(
                np.isnan(model_identifiers[name]),
                np.isnan(static_identifiers[name])))
            self.assertTrue(np.allclose(
                np.nan_to_num(model_identifiers[name]),
                np.nan_to_num(static_identifiers[name])))

    def test_group_identifiers_reject_invalid_shapes(self):
        """Static group identifiers require two three-component bounds."""
        self.assertRaises(
            ValueError, Boundary.group_coordinate_identifiers,
            [0.0, 0.0], [1.0, 1.0, 1.0])

    def test_group_identifiers_reject_nonfinite_bounds(self):
        """Static group identifiers reject nonfinite coordinates."""
        self.assertRaises(
            ValueError, Boundary.group_coordinate_identifiers,
            [0.0, 0.0, 0.0], [1.0, np.inf, 1.0])

    def test_group_identifiers_reject_reversed_bounds(self):
        """Static group identifiers reject minima above maxima."""
        self.assertRaises(
            ValueError, Boundary.group_coordinate_identifiers,
            [0.0, 0.0, 0.0], [1.0, -1.0, 1.0])

    def test_group_identifiers_reject_invalid_tolerance(self):
        """Static group identifiers require a nonnegative tolerance."""
        self.assertRaises(
            ValueError, Boundary.group_coordinate_identifiers,
            [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], -1.0)
        self.assertRaises(
            ValueError, Boundary.group_coordinate_identifiers,
            [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], None)

    def test_constructor_rejects_nonpositive_sampling_density(self):
        """Boundary sampling requires a strictly positive density."""
        fixtures.simple_cubic("DENSITY", length=2.0)
        self.assertRaises(ValueError, Boundary, "DENSITY",
                          sampling_density=0)
        self.assertRaises(ValueError, Boundary, "DENSITY",
                          sampling_density=1.5)
        self.assertRaises(ValueError, Boundary, "DENSITY",
                          sampling_density=True)


class TestBoundaryClassification(GmshTestCase):

    def test_corner_is_classified_into_the_correct_point_group(self):
        """The minimum corner is assigned to POINT_1."""
        fixtures.simple_cubic("CC", length=2.0)
        boundary = Boundary("CC")
        dim, tag = boundary["POINT_1"][0]
        coord = gmsh.model.get_value(dim, tag, [])
        self.assertTrue(np.allclose(coord, [-1, -1, -1]))

    def test_corner_cell_populates_point_groups(self):
        """A corner-connected cubic cell populates all point groups."""
        fixtures.simple_cubic("TC", length=1.0)
        boundary = Boundary("TC")
        self.assertTrue(boundary.periodic)
        self.assertEqual(len(_nonempty(boundary, POINT_GROUPS)), 8)

    def test_occ_padding_is_independent_of_geometry_tolerance(self):
        """OCC bounding-box padding does not reject a valid cube."""
        original_tolerance = gmsh.option.get_number("Geometry.Tolerance")
        try:
            gmsh.option.set_number("Geometry.Tolerance", 1.0e-12)
            fixtures.solid_cube("TIGHT_TOLERANCE", length=2.0)
            boundary = Boundary("TIGHT_TOLERANCE")
            self.assertTrue(boundary.periodic)
        finally:
            gmsh.option.set_number(
                "Geometry.Tolerance", original_tolerance)

    def test_noncorner_cell_populates_edge_groups(self):
        """A non-corner cell populates edge rather than point groups."""
        fixtures.plus_lattice_2d("TP", size=1.0)
        boundary = Boundary("TP")
        self.assertTrue(boundary.periodic)
        self.assertEqual(len(_nonempty(boundary, POINT_GROUPS)), 0)
        self.assertEqual(len(_nonempty(boundary, EDGE_GROUPS[0:4])), 4)

    def test_trimmed_sphere_classifies_and_meshes_periodically(self):
        """A cell with degenerate internal OCC curves remains periodic."""
        fixtures.spherical_inclusion_cube(
            "SPHERE", length=2.0, radius=0.4)
        cell = GmshRVE("SPHERE")
        boundary = cell.boundary

        self.assertTrue(boundary.periodic)
        self.assertEqual(
            sum(len(boundary[name]) for name in POINT_GROUPS), 8)
        self.assertEqual(
            sum(len(boundary[name]) for name in EDGE_GROUPS), 12)
        self.assertEqual(
            sum(len(boundary[name]) for name in FACE_GROUPS), 6)

        cell.add_periodic_mesh_constraints()
        cell.generate_mesh(element_size=0.4)
        cell.build_node_sets()
        self.assertTrue(cell.periodic_mesh)
        for secondary, primary in PERIODIC_PAIRS.items():
            self.assertEqual(len(cell.node_sets[secondary]),
                             len(cell.node_sets[primary]))

    def test_fragmented_exterior_entities_are_classified_and_paired(self):
        """Multiple entities per cuboidal face remain periodic."""
        fixtures.fragmented_cube("FRAGMENTED", length=2.0)
        boundary = Boundary("FRAGMENTED")

        self.assertTrue(boundary.periodic)
        self.assertEqual(len(_nonempty(boundary, FACE_GROUPS)), 6)
        self.assertGreater(max(len(boundary[name]) for name in FACE_GROUPS),
                           1)

    def test_lower_dimensional_entities_remain_in_face_groups(self):
        """Nonconforming beam endpoints on cell faces are retained."""
        fixtures.matrix_with_boundary_fiber("FACE_POINTS", size=1.0)
        boundary = Boundary("FACE_POINTS")

        self.assertTrue(boundary.periodic)
        self.assertTrue(any(dim == 0 for dim, _ in boundary["FACE_1"]))
        self.assertTrue(any(dim == 0 for dim, _ in boundary["FACE_2"]))


class TestBoundaryPeriodicity(GmshTestCase):

    def test_secondary_entities_without_primary_partner_are_rejected(self):
        """A secondary-only group cannot be discarded as periodic."""
        fixtures.square_without_primary_edge("SECONDARY_ONLY", size=2.0)
        boundary = Boundary("SECONDARY_ONLY")

        self.assertFalse(boundary.periodic)
        self.assertEqual(len(boundary["EDGE_4"]), 0)
        self.assertEqual(len(boundary["EDGE_2"]), 1)

    def test_secondary_entity_cannot_match_multiple_primary_entities(self):
        """Periodic entity pairing is bijective."""
        fixtures.square_with_duplicate_primary_edge(
            "DUPLICATE_PRIMARY", size=2.0)
        boundary = Boundary("DUPLICATE_PRIMARY")

        self.assertFalse(boundary.periodic)
        self.assertEqual(len(boundary["EDGE_4"]), 2)
        self.assertEqual(len(boundary["EDGE_2"]), 1)

    def test_empty_samples_do_not_establish_periodicity(self):
        """Entities without valid samples cannot form a periodic pair."""

        class EmptySamplingBoundary(Boundary):

            @staticmethod
            def _generate_sample_points_on_entity(dim, tag, n, tol):
                """Return no samples for positive-dimensional entities."""
                if dim > 0:
                    return np.empty((0, 3))
                return Boundary._generate_sample_points_on_entity(
                    dim, tag, n, tol)

        fixtures.square_lattice("EMPTY_SAMPLES", size=2.0)
        unpaired = Boundary("EMPTY_SAMPLES", check_periodic=False)
        boundary = EmptySamplingBoundary(
            "EMPTY_SAMPLES", initial_data=unpaired.data)

        self.assertFalse(boundary.periodic)

    def test_periodic_pair_entities_are_offset_by_the_lattice_vector(self):
        """Paired entities differ by the expected lattice vector."""
        length = 2.0
        fixtures.simple_cubic("PP", length=length)
        boundary = Boundary("PP")
        self.assertTrue(boundary.periodic)

        # POINT_2 = (xmax, ymin, zmin) is the x-shifted image of POINT_1
        primary_dim, primary_tag = boundary["POINT_1"][0]
        secondary_dim, secondary_tag = boundary["POINT_2"][0]
        primary = np.array(
            gmsh.model.get_value(primary_dim, primary_tag, []))
        secondary = np.array(
            gmsh.model.get_value(secondary_dim, secondary_tag, []))
        self.assertTrue(np.allclose(secondary - primary,
                                    [length, 0.0, 0.0]))


class TestBoundaryDict(GmshTestCase):

    def test_setitem_rejects_invalid_keys_and_uppercases(self):
        """Boundary keys are validated and normalized to uppercase."""
        fixtures.simple_cubic("KV", length=2.0)
        boundary = Boundary("KV")
        self.assertRaises(KeyError, boundary.__setitem__, "NOPE", [])
        boundary["point_1"] = [(0, 123)]           # lowercase is accepted
        self.assertIn((0, 123), boundary["POINT_1"])  # and stored uppercased


class TestGmshRVEValidation(GmshTestCase):

    def test_two_dimensional_rve_must_use_xy_coordinates(self):
        """A two-dimensional RVE in the x-z plane is rejected."""
        fixtures.square_lattice_in_plane(
            "XZ_PLANE", size=2.0, axes=(0, 2))

        with self.assertRaises(RuntimeError) as context:
            GmshRVE("XZ_PLANE")
        self.assertIn("x-y plane", str(context.exception))

    def test_two_dimensional_rve_must_have_zero_z_coordinate(self):
        """A two-dimensional RVE offset from the x-y plane is rejected."""
        fixtures.square_lattice_in_plane(
            "OFFSET_XY_PLANE", size=2.0, axes=(0, 1), offset=1.0)

        with self.assertRaises(RuntimeError) as context:
            GmshRVE("OFFSET_XY_PLANE")
        self.assertIn("z = 0", str(context.exception))

    def test_reserved_boundary_group_name_is_rejected(self):
        """User groups may not reuse reserved boundary names."""
        fixtures.square_lattice("CONFLICT", size=1.0)
        curve = gmsh.model.get_entities(1)[0][1]
        group = gmsh.model.add_physical_group(1, [curve])
        gmsh.model.set_physical_name(1, group, "POINT_1")
        self.assertRaises(RuntimeError, GmshRVE, "CONFLICT")

    def test_one_dimensional_model_is_rejected(self):
        """GmshRVE rejects models with geometric dimension below two."""
        gmsh.model.add("LINE")
        g = gmsh.model.geo
        line = g.add_line(g.add_point(0.0, 0.0, 0.0),
                          g.add_point(1.0, 0.0, 0.0))
        g.synchronize()
        grp = gmsh.model.add_physical_group(1, [line])
        gmsh.model.set_physical_name(1, grp, "STRUTS")
        self.assertRaises(RuntimeError, GmshRVE, "LINE")

    def test_missing_model_is_rejected(self):
        """GmshRVE rejects a model name that is not present in Gmsh."""
        self.assertRaises(ValueError, GmshRVE, "DOES_NOT_EXIST")
