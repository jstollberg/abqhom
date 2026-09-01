"""
Tests for abqhom.utils (Voigt/tensor conversion and the analytical
Voigt/Reuss bounds). These do not need Gmsh or Abaqus.

"""
import unittest

import numpy as np

from abqhom.utils import voigt_to_tensor, tensor_to_voigt, voigt, reuss


class TestVoigtTensor(unittest.TestCase):

    def test_roundtrip_3d_strain(self):
        """voigt -> tensor -> voigt reproduces a 3d strain vector."""
        v = np.array([0.1, -0.2, 0.3, 0.4, -0.5, 0.6])
        back = tensor_to_voigt(voigt_to_tensor(v, "STRAIN"), "STRAIN")
        self.assertTrue(np.allclose(back, v))

    def test_roundtrip_2d_stress(self):
        """voigt -> tensor -> voigt reproduces a 2d stress vector."""
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
        # Two isotropic phases in Voigt (3d) notation
        self.C1 = np.diag([10.0, 10.0, 10.0, 4.0, 4.0, 4.0])
        self.C2 = np.diag([1.0, 1.0, 1.0, 0.4, 0.4, 0.4])
        self.fraction = 0.3

    def test_voigt_is_rule_of_mixtures(self):
        """The Voigt estimate is the volume-weighted average."""
        expected = self.fraction*self.C1 + (1 - self.fraction)*self.C2
        self.assertTrue(np.allclose(voigt(self.fraction, self.C1, self.C2),
                                    expected))

    def test_reuss_below_voigt(self):
        """The Reuss estimate is not stiffer than the Voigt estimate."""
        C_voigt = voigt(self.fraction, self.C1, self.C2)
        C_reuss = reuss(self.fraction, self.C1, self.C2)
        eig_voigt = np.sort(np.linalg.eigvals(C_voigt).real)
        eig_reuss = np.sort(np.linalg.eigvals(C_reuss).real)
        self.assertTrue(np.all(eig_reuss <= eig_voigt + 1e-9))


if __name__ == "__main__":
    unittest.main()
