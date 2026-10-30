"""
Granular tests for the Boundary class and the GmshRVE input validation.
Gmsh only, no Abaqus kernel.

"""
import numpy as np

import gmsh

from abqhom.gmshmodel import GmshRVE, Boundary
from abqhom.constants import POINT_GROUPS, EDGE_GROUPS
from abqhom.tests import fixtures
from abqhom.tests.base import GmshTestCase


def _nonempty(boundary, names):
    """Return the supplied boundary groups that contain entities."""
    result = []
    for name in names:
        if len(boundary[name]) > 0:
            result.append(name)
    return result


class TestBoundaryGeometry(GmshTestCase):

    def test_bounding_box_from_points(self):
        """Boundary points recover the exact cuboidal bounding box."""
        fixtures.simple_cubic("BB", length=2.0)
        box = Boundary("BB").get_bounding_box_from_points()
        self.assertTrue(np.allclose(box, [-1, -1, -1, 1, 1, 1]))

    def test_group_coordinate_identifiers(self):
        """Group identifiers encode fixed and unconstrained coordinates."""
        fixtures.simple_cubic("ID", length=2.0)
        identifiers = Boundary("ID").get_group_coordinate_identifiers()
        # POINT_1 is the min corner, POINT_7 the max corner
        self.assertTrue(np.allclose(identifiers["POINT_1"], [-1, -1, -1]))
        self.assertTrue(np.allclose(identifiers["POINT_7"], [1, 1, 1]))
        # FACE_1 lies at z = zmin and is free in x and y
        face1 = identifiers["FACE_1"]
        self.assertTrue(np.isnan(face1[0]) and np.isnan(face1[1]))
        self.assertAlmostEqual(face1[2], -1.0)

    def test_constructor_rejects_nonpositive_sampling_density(self):
        """Boundary sampling requires a strictly positive density."""
        fixtures.simple_cubic("DENSITY", length=2.0)
        self.assertRaises(ValueError, Boundary, "DENSITY",
                          sampling_density=0)


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

    def test_noncorner_cell_populates_edge_groups(self):
        """A cell that crosses the boundary away from the corners stays
        periodic, with empty POINT groups and populated EDGE groups."""
        fixtures.plus_lattice_2d("TP", size=1.0)
        boundary = Boundary("TP")
        self.assertTrue(boundary.periodic)
        self.assertEqual(len(_nonempty(boundary, POINT_GROUPS)), 0)
        self.assertEqual(len(_nonempty(boundary, EDGE_GROUPS[0:4])), 4)


class TestBoundaryPeriodicity(GmshTestCase):

    def test_periodic_pair_entities_are_offset_by_the_lattice_vector(self):
        """A primary entity and its paired secondary entity (same list index)
        differ by exactly the periodicity offset."""
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
