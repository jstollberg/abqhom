"""Tests for the explicit-path Gmsh loading contract of the test runner."""

import os
import shutil
import tempfile
import unittest

import gmsh

from abqhom.tests import _load_gmsh


class TestGmshLoading(unittest.TestCase):

    def test_missing_path_is_rejected(self):
        """The runner requires the caller to prescribe a Gmsh path."""
        self.assertRaises(ValueError, _load_gmsh, None)

    def test_directory_without_module_is_rejected(self):
        """The prescribed directory must contain the Gmsh Python module."""
        missing = os.path.join(os.getcwd(), "missing_gmsh_directory")
        self.assertRaises(IOError, _load_gmsh, missing)

    def test_loaded_module_from_prescribed_path_is_accepted(self):
        """An already loaded Gmsh module is accepted only at its own path."""
        gmsh_path = os.path.dirname(gmsh.__file__)
        self.assertIs(_load_gmsh(gmsh_path), gmsh)

    def test_loaded_module_from_other_path_is_rejected(self):
        """A loaded Gmsh module cannot be replaced by another installation."""
        other = tempfile.mkdtemp(prefix="other_gmsh_")
        try:
            with open(os.path.join(other, "gmsh.py"), "w"):
                pass
            self.assertRaises(RuntimeError, _load_gmsh, other)
        finally:
            shutil.rmtree(other)
