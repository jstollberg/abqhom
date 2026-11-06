# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-or-later

"""
End-to-end homogenization tests for AbaqusRVE. These must be run in the
Abaqus Python interpreter (they build a small RVE, run a linear analysis
and check the effective stiffness). Meshes are kept tiny for the
Learning Edition node limit.

"""
import numpy as np

from abqhom.utils import voigt_to_tensor, tensor_to_voigt
from abqhom.gmshmodel import GmshRVE
from abqhom.abqmodel import AbaqusRVE, LoadHistory
from abqhom.tests import fixtures
from abqhom.tests.base import AbaqusTestCase

STRAIN = 1.0e-3


def _stiffness_column(build, assign_sections, comp, element_size,
                      custom_types=None, bc="periodic"):
    """Return one effective-stiffness column from a unit-strain case."""
    name = build()
    rve = GmshRVE(name)
    rve.add_periodic_mesh_constraints()
    rve.generate_mesh(element_size=element_size)

    unit_cell = AbaqusRVE(rve, custom_types=custom_types)
    assign_sections(unit_cell)

    def strain_history(time):
        """Return the selected ramped macroscopic strain component."""
        strain = np.zeros(6)
        strain[comp] = STRAIN
        return time*voigt_to_tensor(strain, "STRAIN")

    unit_cell.model.StaticStep(name="LOAD", previous="Initial")
    history = LoadHistory(np.array([0.0, 1.0]), strain_history)
    if bc == "affine":
        unit_cell.apply_affine_displacement_bc("LOAD", history)
    else:
        unit_cell.apply_periodic_bc("LOAD", history)
    unit_cell.run_analysis(field_outputs=["U"])
    stress = tensor_to_voigt(unit_cell.average_stress("LOAD"), "STRESS")
    return stress/STRAIN


def _elastic(unit_cell, name, E, nu):
    """Create an isotropic linear-elastic Abaqus material."""
    material = unit_cell.model.Material(name=name)
    material.Elastic(table=[[E, nu]])
    return material


class TestBeamLattices(AbaqusTestCase):
    """Test beam-lattice stiffness and incremental power."""

    def test_simple_cubic_axial_modulus(self):
        """The simple-cubic axial modulus matches its analytical value."""
        E, nu, r, L = 210.0, 0.3, 0.05, 1.0

        def assign(unit_cell):
            """Assign the linear-elastic beam section."""
            _elastic(unit_cell, "M", E, nu)
            unit_cell.assign_beam_section("STRUTS", "M", r)

        column = _stiffness_column(
            lambda: fixtures.simple_cubic("SC", length=L), assign, 0, L)
        analytical = E*np.pi*r**2/L**2
        self.assertAlmostEqual(column[0], analytical,
                               delta=1.0e-3*analytical)
        self.assertLess(abs(column[1]), 1.0e-3*analytical)

    def test_boundary_strut_ownership(self):
        """Reduced square-edge struts recover the owned axial stiffness."""
        E, nu, r = 1.0, 0.3, 0.05
        area = np.pi*r*r

        def assign(unit_cell):
            """Assign the linear-elastic beam section."""
            _elastic(unit_cell, "M", E, nu)
            unit_cell.assign_beam_section("STRUTS", "M", r)

        column = _stiffness_column(
            lambda: fixtures.square_lattice("SQ", size=1.0), assign, 0, 0.5)
        self.assertAlmostEqual(column[0], E*area, delta=1.0e-4*E*area)

    def test_2d_virtual_power_matches_macroscopic_power(self):
        """Microscopic and macroscopic power agree for a planar beam cell."""
        E, nu, r = 1.0, 0.3, 0.05
        fixtures.square_lattice("VP2", size=1.0)
        rve = GmshRVE("VP2")
        rve.add_periodic_mesh_constraints()
        rve.generate_mesh(element_size=1.0)

        unit = AbaqusRVE(rve)
        _elastic(unit, "M", E, nu)
        unit.assign_beam_section("STRUTS", "M", r)
        unit.model.StaticStep(name="LOAD", previous="Initial")

        def strain_history(time):
            """Return the ramped axial strain history."""
            strain = np.zeros(6)
            strain[0] = time*STRAIN
            return voigt_to_tensor(strain, "STRAIN")

        history = LoadHistory(np.array([0.0, 1.0]), strain_history)
        unit.apply_periodic_bc("LOAD", history)
        unit.run_analysis()

        stress = unit.average_stress("LOAD")
        macroscopic_power = float(np.sum(
            stress*strain_history(1.0)[:2, :2]))
        microscopic_power = unit.average_power_increment("LOAD")

        self.assertAlmostEqual(microscopic_power, macroscopic_power,
                               delta=1.0e-6*abs(macroscopic_power))

    def test_3d_rotational_power_matches_macroscopic_power(self):
        """Beam moments and finite rotation increments enter 3D power."""
        E, nu, r = 10.0, 0.3, 0.05
        fixtures.simple_cubic("VP3", length=1.0)
        rve = GmshRVE("VP3")
        rve.add_periodic_mesh_constraints()
        rve.generate_mesh(element_size=0.5)

        unit = AbaqusRVE(rve)
        _elastic(unit, "M", E, nu)
        unit.assign_beam_section("STRUTS", "M", r)
        unit.model.StaticStep(name="LOAD", previous="Initial")

        def strain_history(time):
            """Return the ramped engineering-shear history."""
            strain = np.zeros(6)
            strain[5] = time*STRAIN
            return voigt_to_tensor(strain, "STRAIN")

        history = LoadHistory(np.array([0.0, 1.0]), strain_history)
        unit.apply_periodic_bc("LOAD", history)
        unit.run_analysis()

        stress = unit.average_stress("LOAD")
        macroscopic_power = float(np.sum(stress*strain_history(1.0)))
        microscopic_power = unit.average_power_increment("LOAD")
        self.assertAlmostEqual(microscopic_power, macroscopic_power,
                               delta=1.0e-6*abs(macroscopic_power))

    def test_bcc_axis_permutation_and_major_symmetry(self):
        """The BCC response satisfies C11=C22 and major symmetry C12=C21."""
        E, nu, r, L = 210.0, 0.3, 0.05, 1.0

        def assign(unit_cell):
            """Assign identical sections to edge and diagonal beams."""
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

    def test_one_model_runs_several_load_cases(self):
        """One AbaqusRVE supports consecutive periodic load cases."""
        E, nu, r, L = 210.0, 0.3, 0.05, 1.0
        fixtures.simple_cubic("REUSE", length=L)
        rve = GmshRVE("REUSE")
        rve.add_periodic_mesh_constraints()
        rve.generate_mesh(element_size=L)
        unit = AbaqusRVE(rve)
        _elastic(unit, "M", E, nu)
        unit.assign_beam_section("STRUTS", "M", r)
        unit.model.StaticStep(name="LOAD", previous="Initial")

        columns = []
        for comp in (0, 1):
            def strain_history(time, c=comp):
                """Return the current ramped normal-strain history."""
                strain = np.zeros(6)
                strain[c] = STRAIN
                return time*voigt_to_tensor(strain, "STRAIN")

            unit.apply_periodic_bc("LOAD", LoadHistory(np.array([0.0, 1.0]),
                                                       strain_history))
            unit.run_analysis(field_outputs=["U"])
            columns.append(tensor_to_voigt(unit.average_stress("LOAD"),
                                           "STRESS")
                           / STRAIN)

        self.assertEqual(unit.job_name, "JOB_REUSE_2")
        analytical = E*np.pi*r**2/L**2
        self.assertAlmostEqual(columns[0][0], analytical,
                               delta=1.0e-3*analytical)
        self.assertAlmostEqual(columns[0][0], columns[1][1],
                               delta=1.0e-6*columns[0][0])


class TestShells(AbaqusTestCase):
    """Test homogenization of periodic shell structures."""

    def test_parallel_plates_membrane_modulus(self):
        """Periodic shell plates recover their analytical membrane moduli."""
        E, nu, t, size, spacing = 210.0, 0.3, 0.01, 1.0, 1.0

        def assign(unit_cell):
            """Assign the linear-elastic shell section."""
            _elastic(unit_cell, "M", E, nu)
            unit_cell.assign_shell_section("SHELLS", "M", t)

        column = _stiffness_column(
            lambda: fixtures.parallel_plates_3d("SH", size=size,
                                                spacing=spacing),
            assign, 0, 0.5*size)
        C11 = E*t/((1.0 - nu**2)*spacing)
        self.assertAlmostEqual(column[0], C11, delta=1.0e-4*C11)
        self.assertAlmostEqual(column[1], nu*C11, delta=1.0e-4*C11)
        self.assertLess(abs(column[2]), 1.0e-4*C11)


class TestContinuum(AbaqusTestCase):
    """Test periodic and affine homogenization of a continuum cell."""

    def test_homogeneous_cube_recovers_material(self):
        """A homogeneous continuum RVE homogenizes to the material itself."""
        E, nu, L = 100.0, 0.3, 1.0

        def assign(unit_cell):
            """Assign the homogeneous solid section."""
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

    def test_affine_bc_recovers_material_for_homogeneous_cube(self):
        """Affine BCs recover a homogeneous continuum's material stiffness."""
        E, nu, L = 100.0, 0.3, 1.0

        def assign(unit_cell):
            """Assign the homogeneous solid section."""
            _elastic(unit_cell, "M", E, nu)
            unit_cell.assign_solid_section("SOLID", "M")

        factor = E/((1 + nu)*(1 - 2*nu))
        C11, C12 = factor*(1 - nu), factor*nu
        C44 = E/(2*(1 + nu))

        cxx = _stiffness_column(
            lambda: fixtures.solid_cube("CBA", length=L), assign, 0, 0.5*L,
            bc="affine")
        cyz = _stiffness_column(
            lambda: fixtures.solid_cube("CBA", length=L), assign, 3, 0.5*L,
            bc="affine")
        self.assertAlmostEqual(cxx[0], C11, delta=1.0e-3*C11)
        self.assertAlmostEqual(cxx[1], C12, delta=1.0e-3*C11)
        self.assertAlmostEqual(cyz[3], C44, delta=1.0e-3*C44)


class TestMixedDimensionalCell(AbaqusTestCase):
    """Test reduction and power equivalence of a mixed cell."""

    def test_reduction_and_incremental_power_equivalence(self):
        """A reduced continuum--shell--beam cell satisfies Hill--Mandel."""
        E, nu = 100.0, 0.3
        fixtures.mixed_dimensional_cube("MIXED", size=1.0,
                                        nodes_per_edge=3)
        rve = GmshRVE("MIXED")
        rve.add_periodic_mesh_constraints()
        rve.generate_mesh()
        rve.build_element_sets()

        expected_elements = (len(rve.element_sets["ALL_ELEMENTS"])
                             - len(rve.element_sets["SECONDARY_ELEMENTS"]))
        self.assertGreater(len(rve.element_sets["SECONDARY_ELEMENTS"]), 0)

        unit = AbaqusRVE(
            rve, custom_types={1: "B31", 3: "S4", 5: "C3D8"})
        self.assertEqual(len(unit.part.elements), expected_elements)
        for name in ("CONTINUUM", "SHELLS", "BEAMS"):
            self.assertGreater(
                len(unit.part.sets["ELSET_" + name].elements), 0)

        _elastic(unit, "M", E, nu)
        unit.assign_solid_section("CONTINUUM", "M")
        unit.assign_shell_section("SHELLS", "M", 0.02)
        unit.assign_beam_section("BEAMS", "M", 0.04)
        unit.model.StaticStep(name="LOAD", previous="Initial")

        def strain_history(time):
            """Apply a small normal macroscopic displacement gradient."""
            gradient = np.zeros((3, 3))
            gradient[0, 0] = time*STRAIN
            return gradient

        history = LoadHistory(np.array([0.0, 1.0]), strain_history)
        unit.apply_periodic_bc("LOAD", history)
        unit.run_analysis()

        stress = unit.average_stress("LOAD")
        macroscopic_power = float(np.sum(stress*strain_history(1.0)))
        microscopic_power = unit.average_power_increment("LOAD")
        scale = max(abs(macroscopic_power), abs(microscopic_power), 1.0e-14)

        self.assertTrue(np.all(np.isfinite(stress)))
        self.assertGreater(stress[0, 0], 0.0)
        self.assertLess(abs(microscopic_power - macroscopic_power)/scale,
                        1.0e-6)


class TestEmbeddedFibers(AbaqusTestCase):
    """Test periodic homogenization with embedded beam fibers."""

    def test_boundary_fiber_excluded_from_translation_bc(self):
        """Embedded boundary fibers are excluded from periodic translations."""
        E_m, nu_m, E_f, nu_f = 50.0, 0.3, 210.0, 0.2

        fixtures.matrix_with_boundary_fiber("EMB", size=1.0)
        rve = GmshRVE("EMB")
        rve.add_periodic_mesh_constraints()
        rve.generate_mesh()

        unit = AbaqusRVE(rve)
        _elastic(unit, "M", E_m, nu_m)
        _elastic(unit, "F", E_f, nu_f)
        unit.assign_solid_section("MATRIX", "M")
        unit.assign_beam_section("FIBER", "F", 0.1)
        unit.model.EmbeddedRegion("COUPLING",
                                  unit.instance.sets["ELSET_FIBER"],
                                  unit.instance.sets["ELSET_MATRIX"])
        unit.model.StaticStep(name="LOAD", previous="Initial")

        def strain_history(time):
            """Return the ramped axial strain applied to the fiber."""
            strain = np.zeros(6)
            strain[2] = STRAIN
            return time*voigt_to_tensor(strain, "STRAIN")

        unit.apply_periodic_bc("LOAD", LoadHistory(np.array([0.0, 1.0]),
                                                   strain_history))
        unit.run_analysis(field_outputs=["U"])
        stress = tensor_to_voigt(unit.average_stress("LOAD"), "STRESS")

        self.assertTrue(np.all(np.isfinite(stress)))
        constrained_matrix_modulus = (
            E_m*(1.0 - nu_m)/((1.0 + nu_m)*(1.0 - 2.0*nu_m)))
        self.assertGreater(stress[2]/STRAIN, constrained_matrix_modulus)
