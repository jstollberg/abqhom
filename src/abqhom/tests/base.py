# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-only

"""Shared test-case bases for Gmsh and Abaqus tests."""

import unittest

import gmsh


class GmshTestCase(unittest.TestCase):
    """Initialize Gmsh and remove all models after every test."""

    @classmethod
    def setUpClass(cls):
        """Initialize the prescribed Gmsh module once per test class."""
        if not gmsh.is_initialized():
            gmsh.initialize()

    def tearDown(self):
        """Remove every Gmsh model created by the completed test."""
        for name in list(gmsh.model.list()):
            gmsh.model.set_current(name)
            gmsh.model.remove()


class AbaqusTestCase(GmshTestCase):
    """Close output databases and remove Gmsh models after every test."""

    def tearDown(self):
        """Release Abaqus resources before removing the Gmsh fixtures."""
        try:
            from abaqus import session
            for name in list(session.odbs.keys()):
                session.odbs[name].close()
        finally:
            super(AbaqusTestCase, self).tearDown()
