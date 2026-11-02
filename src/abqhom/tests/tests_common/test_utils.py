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

    def test_roundtrip_3d_strain(self):
        """Voigt to tensor to Voigt reproduces a 3D strain vector."""
        v = np.array([0.1, -0.2, 0.3, 0.4, -0.5, 0.6])
        back = tensor_to_voigt(voigt_to_tensor(v, "STRAIN"), "STRAIN")
        self.assertTrue(np.allclose(back, v))

    def test_roundtrip_2d_stress(self):
        """Voigt to tensor to Voigt reproduces a 2D stress vector."""
        v = np.array([1.0, -2.0, 0.5])
        back = tensor_to_voigt(voigt_to_tensor(v, "STRESS"), "STRESS",
                               check_symmetry=False)
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


class TestBounds(unittest.TestCase):

    def setUp(self):
        # Two diagonal positive-definite stiffness matrices
        self.C1 = np.diag([10.0, 10.0, 10.0, 4.0, 4.0, 4.0])
        self.C2 = np.diag([1.0, 1.0, 1.0, 0.4, 0.4, 0.4])
        self.fraction = 0.3

    def test_voigt_is_rule_of_mixtures(self):
        """The Voigt estimate is the volume-weighted average."""
        expected = self.fraction*self.C1 + (1 - self.fraction)*self.C2
        self.assertTrue(np.allclose(voigt(self.fraction, self.C1, self.C2),
                                    expected))

    def test_reuss_matches_harmonic_mixture_and_is_below_voigt(self):
        """The Reuss estimate is the harmonic mixture and a lower bound."""
        C_voigt = voigt(self.fraction, self.C1, self.C2)
        C_reuss = reuss(self.fraction, self.C1, self.C2)
        expected = np.linalg.inv(
            self.fraction*np.linalg.inv(self.C1)
            + (1.0 - self.fraction)*np.linalg.inv(self.C2))
        self.assertTrue(np.allclose(C_reuss, expected))
        difference = 0.5*((C_voigt - C_reuss)
                          + (C_voigt - C_reuss).T)
        self.assertGreaterEqual(np.min(np.linalg.eigvalsh(difference)),
                                -1.0e-12)


class TestHillMandelResidual(unittest.TestCase):

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

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="abqhom_test_")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_csv_roundtrip(self):
        """export_to_csv_file then read_csv_file recovers the array."""
        data = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        base = os.path.join(self.tmp, "stiffness")
        export_to_csv_file(data, base)
        back = read_csv_file(base + ".csv")
        self.assertTrue(np.allclose(back, data))
