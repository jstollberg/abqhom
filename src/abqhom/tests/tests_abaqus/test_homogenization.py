"""
End-to-end homogenization tests for AbaqusRVE. These must be run in the
Abaqus Python interpreter (they build a small RVE, run a linear analysis
and check the effective stiffness). Meshes are kept tiny for the
Learning Edition node limit.

"""
import unittest

import numpy as np

import gmsh

from abqhom.utils import voigt_to_tensor, tensor_to_voigt
from abqhom.gmshmodel import GmshRVE
from abqhom.abqmodel import AbaqusRVE, LoadHistory
from abqhom.tests import fixtures

STRAIN = 1.0e-3


def _stiffness_column(build, assign_sections, comp, element_size,
                      custom_types=None):
    """Run one unit-strain load case and return the corresponding column
    of the effective stiffness in Voigt notation."""
    name = build()
    rve = GmshRVE(name)
    rve.add_periodic_mesh_constraints()
    rve.generate_mesh(element_size=element_size)

    unit_cell = AbaqusRVE(rve, custom_types or {})
    assign_sections(unit_cell)

    def strain_history(time):
        strain = np.zeros(6)
        strain[comp] = STRAIN
        return time*voigt_to_tensor(strain, "STRAIN")

    unit_cell.model.StaticStep(name="LOAD", previous="Initial")
    unit_cell.apply_periodic_BC("LOAD", LoadHistory(np.array([0.0, 1.0]),
                                                    strain_history))
    unit_cell.run_analysis(field_outputs=["U"])
    stress = tensor_to_voigt(unit_cell.average_stress("LOAD"), "STRESS",
                             check_symmetry=False)

    gmsh.model.set_current(name)
    gmsh.model.remove()
    return stress/STRAIN


def _elastic(unit_cell, name, E, nu):
    material = unit_cell.model.Material(name=name)
    material.Elastic(table=[[E, nu]])
    return material


class TestBeamLattices(unittest.TestCase):

    def test_simple_cubic_axial_modulus(self):
        """C11 of a simple-cubic beam lattice matches E*pi*r**2/L**2 and
        the transverse coupling C12 vanishes."""
        E, nu, r, L = 210.0, 0.3, 0.05, 1.0

        def assign(unit_cell):
            _elastic(unit_cell, "M", E, nu)
            unit_cell.assign_beam_section("STRUTS", "M", r)

        column = _stiffness_column(
            lambda: fixtures.simple_cubic("SC", length=L), assign, 0, L)
        analytical = E*np.pi*r**2/L**2
        self.assertAlmostEqual(column[0], analytical,
                               delta=1.0e-3*analytical)
        self.assertLess(abs(column[1]), 1.0e-3*analytical)

    def test_boundary_strut_sharing_weight(self):
        """The four edge struts of a 2d square cell are shared by two cells
        each, so C11 = E*A (not 2*E*A) thanks to the sharing weights."""
        E, nu, r = 1.0, 0.3, 0.05
        area = np.pi*r*r

        def assign(unit_cell):
            _elastic(unit_cell, "M", E, nu)
            unit_cell.assign_beam_section("STRUTS", "M", r)

        column = _stiffness_column(
            lambda: fixtures.square_lattice("SQ", size=1.0), assign, 0, 1.0)
        self.assertAlmostEqual(column[0], E*area, delta=1.0e-4*E*area)

    def test_bcc_cubic_symmetry(self):
        """The BCC lattice has cubic symmetry (C11 = C22) and the effective
        stiffness is minor-symmetric (C12 = C21)."""
        E, nu, r, L = 210.0, 0.3, 0.05, 1.0

        def assign(unit_cell):
            _elastic(unit_cell, "M", E, nu)
            unit_cell.assign_beam_section("CUBIC", "M", r)
            unit_cell.assign_beam_section("BCC", "M", r)

        cxx = _stiffness_column(
            lambda: fixtures.bcc("BC", length=L), assign, 0, L)
        cyy = _stiffness_column(
            lambda: fixtures.bcc("BC", length=L), assign, 1, L)
        self.assertAlmostEqual(cxx[0], cyy[1], delta=1.0e-6*abs(cxx[0]))
        self.assertAlmostEqual(cxx[1], cyy[0],
                               delta=1.0e-6*abs(cxx[0]) + 1.0e-9)


class TestContinuum(unittest.TestCase):

    def test_homogeneous_cube_recovers_material(self):
        """A homogeneous continuum RVE homogenizes to the material itself."""
        E, nu, L = 100.0, 0.3, 1.0

        def assign(unit_cell):
            _elastic(unit_cell, "M", E, nu)
            unit_cell.assign_solid_section("SOLID", "M")

        factor = E/((1 + nu)*(1 - 2*nu))
        C11, C12 = factor*(1 - nu), factor*nu
        C44 = E/(2*(1 + nu))

        cxx = _stiffness_column(
            lambda: fixtures.solid_cube("CB", length=L), assign, 0, 0.5*L)
        cyz = _stiffness_column(
            lambda: fixtures.solid_cube("CB", length=L), assign, 3, 0.5*L)
        self.assertAlmostEqual(cxx[0], C11, delta=1.0e-3*C11)
        self.assertAlmostEqual(cxx[1], C12, delta=1.0e-3*C11)
        self.assertAlmostEqual(cyz[3], C44, delta=1.0e-3*C44)


if __name__ == "__main__":
    unittest.main()
