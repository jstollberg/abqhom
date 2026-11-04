# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Tests for abqhom.utils (Voigt/tensor conversion and the analytical
Voigt/Reuss bounds). These do not need Gmsh or Abaqus.

"""
import os
import shutil
import tempfile
import unittest

import numpy as np

from abqhom.utils import (voigt_to_tensor, tensor_to_voigt, voigt, reuss,
                          hill_mandel_residual, export_to_csv_file,
                          read_csv_file)


class TestVoigtTensor(unittest.TestCase):
    """Test conversion between tensor and Voigt representations."""

    def test_roundtrip_3d_strain(self):
        """Voigt to tensor to Voigt reproduces a 3D strain vector."""
        v = np.array([0.1, -0.2, 0.3, 0.4, -0.5, 0.6])
        back = tensor_to_voigt(voigt_to_tensor(v, "STRAIN"), "STRAIN")
        self.assertTrue(np.allclose(back, v))

    def test_roundtrip_2d_stress(self):
        """Voigt to tensor to Voigt reproduces a 2D stress vector."""
        v = np.array([1.0, -2.0, 0.5])
        back = tensor_to_voigt(voigt_to_tensor(v, "STRESS"), "STRESS")
        self.assertTrue(np.allclose(back, v))

    def test_strain_shear_halving(self):
        """The engineering shear strain becomes the tensor component/2."""
        tensor = voigt_to_tensor(np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.2]),
                                 "STRAIN")
        self.assertAlmostEqual(tensor[0, 1], 0.1)

    def test_tensor_symmetry_check(self):
        """tensor_to_voigt rejects a non-symmetric stress tensor."""
        unsym = np.array([[1.0, 2.0, 0.0],
                          [0.0, 1.0, 0.0],
                          [0.0, 0.0, 1.0]])
        self.assertRaises(ValueError, tensor_to_voigt, unsym, "STRESS")

    def test_tensor_symmetry_check_can_be_disabled(self):
        """A caller may explicitly convert a nonsymmetric tensor."""
        tensor = np.array([[1.0, 2.0], [3.0, 4.0]])
        result = tensor_to_voigt(
            tensor, "STRESS", check_symmetry=False)
        self.assertTrue(np.array_equal(result, [1.0, 4.0, 2.0]))

    def test_invalid_quantities_are_rejected(self):
        """Both conversion directions reject unknown quantities."""
        self.assertRaises(
            ValueError, voigt_to_tensor, np.zeros(6), "UNKNOWN")
        self.assertRaises(
            ValueError, tensor_to_voigt, np.zeros((3, 3)), "UNKNOWN")

    def test_invalid_shapes_are_rejected(self):
        """Both conversion directions reject unsupported input shapes."""
        self.assertRaises(
            ValueError, voigt_to_tensor, np.zeros(4), "STRESS")
        self.assertRaises(
            ValueError, tensor_to_voigt, np.zeros((4, 4)), "STRESS")


class TestBounds(unittest.TestCase):
    """Test generalized Voigt and Reuss stiffness estimates."""

    def setUp(self):
        """Create three diagonal positive-definite stiffness matrices."""
        self.stiffnesses = np.array([
            np.diag([10.0, 10.0, 10.0, 4.0, 4.0, 4.0]),
            np.diag([4.0, 4.0, 4.0, 1.6, 1.6, 1.6]),
            np.diag([1.0, 1.0, 1.0, 0.4, 0.4, 0.4]),
        ])
        self.fractions = np.array([0.2, 0.3, 0.5])

    def test_voigt_is_rule_of_mixtures(self):
        """The Voigt estimate is the volume-weighted average."""
        expected = sum(fraction*stiffness for fraction, stiffness in zip(
            self.fractions, self.stiffnesses))
        self.assertTrue(np.allclose(
            voigt(self.fractions, self.stiffnesses), expected))

    def test_reuss_matches_harmonic_mixture_and_is_below_voigt(self):
        """The Reuss estimate is the harmonic mixture and a lower bound."""
        C_voigt = voigt(self.fractions, self.stiffnesses)
        C_reuss = reuss(self.fractions, self.stiffnesses)
        expected_compliance = sum(
            fraction*np.linalg.inv(stiffness)
            for fraction, stiffness in zip(
                self.fractions, self.stiffnesses))
        expected = np.linalg.inv(expected_compliance)
        self.assertTrue(np.allclose(C_reuss, expected))
        difference = 0.5*((C_voigt - C_reuss)
                          + (C_voigt - C_reuss).T)
        self.assertGreaterEqual(np.min(np.linalg.eigvalsh(difference)),
                                -1.0e-12)

    def test_invalid_volume_fractions_are_rejected(self):
        """Fractions must be finite, nonnegative, and normalized."""
        self.assertRaises(
            ValueError, voigt, [0.2, 0.3], self.stiffnesses[:2])
        self.assertRaises(
            ValueError, voigt, [0.5, -0.2, 0.7], self.stiffnesses)
        self.assertRaises(
            ValueError, reuss, [0.5, np.inf, 0.5], self.stiffnesses)

    def test_stiffness_data_must_match_constituents(self):
        """Each constituent requires one finite square stiffness matrix."""
        self.assertRaises(
            ValueError, voigt, self.fractions, self.stiffnesses[:2])
        self.assertRaises(
            ValueError, reuss, self.fractions, np.ones((3, 2, 3)))
        invalid = self.stiffnesses.copy()
        invalid[0, 0, 0] = np.nan
        self.assertRaises(ValueError, voigt, self.fractions, invalid)


class TestHillMandelResidual(unittest.TestCase):
    """Test normalization of the discrete Hill--Mandel residual."""

    def test_equal_powers_give_zero(self):
        """Equal incremental-work densities satisfy the condition."""
        residual = hill_mandel_residual(2.0, 2.0, 1.0e-12)
        self.assertEqual(residual, 0.0)

    def test_larger_power_sets_normalization(self):
        """The larger absolute work density determines the normalization."""
        residual = hill_mandel_residual(2.0, 1.0, 0.1)
        self.assertAlmostEqual(residual, 0.5)

    def test_tolerance_regularizes_small_powers(self):
        """The supplied work-density scale regularizes near-zero states."""
        residual = hill_mandel_residual(1.0e-9, 0.0, 1.0e-6)
        self.assertAlmostEqual(residual, 1.0e-3)

    def test_array_inputs_are_evaluated_elementwise(self):
        """Incremental-work histories may be passed as arrays."""
        residual = hill_mandel_residual(np.array([2.0, 1.0e-9]),
                                        np.array([1.0, 0.0]), 1.0e-6)
        self.assertTrue(np.allclose(residual, np.array([0.5, 1.0e-3])))

    def test_nonpositive_tolerance_is_rejected(self):
        """A nonpositive incremental-work scale is not admissible."""
        self.assertRaises(ValueError, hill_mandel_residual, 1.0, 1.0, 0.0)


class TestFileIO(unittest.TestCase):
    """Test CSV export and import utilities."""

    def setUp(self):
        """Create an isolated output directory."""
        self.tmp = tempfile.mkdtemp(prefix="abqhom_test_")

    def tearDown(self):
        """Remove the isolated output directory."""
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_csv_roundtrip(self):
        """export_to_csv_file then read_csv_file recovers the array."""
        data = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        path = os.path.join(self.tmp, "stiffness.txt")
        csv_path = os.path.join(self.tmp, "stiffness.csv")
        export_to_csv_file(data, path)
        back = read_csv_file(csv_path)
        self.assertTrue(np.allclose(back, data))
        self.assertTrue(os.path.isfile(csv_path))
        self.assertFalse(os.path.exists(path))

    def test_csv_roundtrip_with_header(self):
        """CSV headers and custom formats survive a numeric roundtrip."""
        data = np.array([[1.0, 2.5], [3.0, 4.5]])
        path = os.path.join(self.tmp, "values.csv")
        export_to_csv_file(
            data, path, header="first,second", fmt=["%d", "%.1f"])

        back, header = read_csv_file(path, has_header=True)

        self.assertTrue(np.array_equal(back, data))
        self.assertEqual(header, ("first", "second"))

    def test_empty_csv_header_is_rejected(self):
        """Reading a header from an empty CSV file raises a clear error."""
        path = os.path.join(self.tmp, "empty.csv")
        with open(path, "w"):
            pass

        self.assertRaises(
            ValueError, read_csv_file, path, has_header=True)
