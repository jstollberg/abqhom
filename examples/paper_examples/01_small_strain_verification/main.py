# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-or-later

"""Verify small-strain homogenization of a mixed-dimensional unit cell.

The cell combines continuum, shell, and beam constituents. Six independent
strain cases determine its effective stiffness and verify the analytical
normal response, cubic symmetry, the Hill--Mandel condition, and the
equivalence of the reduced and ownership-weighted complete cells.
"""

import os
import sys

import gmsh
import numpy as np

from abqhom.abqmodel import AbaqusRVE, LoadHistory
from abqhom.gmshmodel import GmshRVE
from abqhom.utils import (export_to_csv_file, hill_mandel_residual,
                          tensor_to_voigt, voigt_to_tensor)

here = os.getcwd()


def create_geometry(model_name, edge_length, nodes_per_edge):
    """Create the conforming continuum--shell--beam unit cell."""
    def add_physical_group(dimension, entities, name):
        """Add and name one physical group of the Gmsh model."""
        group = gmsh.model.add_physical_group(dimension, entities)
        gmsh.model.set_physical_name(dimension, group, name)

    if not gmsh.is_initialized():
        gmsh.initialize()
    if model_name in gmsh.model.list():
        gmsh.model.set_current(model_name)
        gmsh.model.remove()
    gmsh.model.add(model_name)

    gmsh.option.set_number("Geometry.Tolerance", 1.0e-8*edge_length)

    h = 0.5*edge_length
    geometry = gmsh.model.geo

    # Four corners at y=-h, followed by the corresponding corners at y=+h
    p1 = geometry.add_point(-h, -h, -h)
    p2 = geometry.add_point(-h, -h, h)
    p3 = geometry.add_point(h, -h, h)
    p4 = geometry.add_point(h, -h, -h)
    p5 = geometry.add_point(-h, h, -h)
    p6 = geometry.add_point(-h, h, h)
    p7 = geometry.add_point(h, h, h)
    p8 = geometry.add_point(h, h, -h)

    lines = [
        geometry.add_line(p1, p2),
        geometry.add_line(p2, p3),
        geometry.add_line(p3, p4),
        geometry.add_line(p4, p1),
        geometry.add_line(p5, p6),
        geometry.add_line(p6, p7),
        geometry.add_line(p7, p8),
        geometry.add_line(p8, p5),
        geometry.add_line(p1, p5),
        geometry.add_line(p2, p6),
        geometry.add_line(p3, p7),
        geometry.add_line(p4, p8),
    ]

    # Surface order: y-, y+, x-, z+, x+, z-
    surface_connectivity = [
        [lines[0], lines[1], lines[2], lines[3]],
        [lines[4], lines[5], lines[6], lines[7]],
        [lines[0], lines[9], lines[4], lines[8]],
        [lines[1], lines[10], lines[5], lines[9]],
        [lines[2], lines[11], lines[6], lines[10]],
        [lines[3], lines[8], lines[7], lines[11]],
    ]
    surfaces = []
    for connectivity in surface_connectivity:
        loop = geometry.add_curve_loop(connectivity, reorient=True)
        surfaces.append(geometry.add_plane_surface([loop]))

    surface_loop = geometry.add_surface_loop(surfaces)
    volume = geometry.add_volume([surface_loop])

    # The transfinite discretization gives all constituents common nodes
    for line in lines:
        geometry.mesh.set_transfinite_curve(line, nodes_per_edge)
    for surface in surfaces:
        geometry.mesh.set_transfinite_surface(surface)
        geometry.mesh.set_recombine(2, surface)
    geometry.mesh.set_transfinite_volume(volume)
    geometry.synchronize()

    add_physical_group(3, [volume], "CONTINUUM")
    add_physical_group(2, [surfaces[2], surfaces[4]], "SHELL_X")
    add_physical_group(2, [surfaces[0], surfaces[1]], "SHELL_Y")
    add_physical_group(2, [surfaces[3], surfaces[5]], "SHELL_Z")
    add_physical_group(1, [lines[1], lines[3], lines[5], lines[7]],
                       "BEAM_X")
    add_physical_group(1, [lines[8], lines[9], lines[10], lines[11]],
                       "BEAM_Y")
    add_physical_group(1, [lines[0], lines[2], lines[4], lines[6]],
                       "BEAM_Z")


def analytical_normal_stiffness(edge_length, continuum, shell, beam,
                                shell_thickness, beam_radius):
    """Return the analytical 3x3 normal-stiffness block."""
    E_c, nu_c = continuum
    E_s, nu_s = shell
    E_b, _ = beam

    mu_c = E_c/(2.0*(1.0 + nu_c))
    lambda_c = E_c*nu_c/((1.0 + nu_c)*(1.0 - 2.0*nu_c))
    shell_normal = E_s/(1.0 - nu_s**2)
    shell_coupling = nu_s*shell_normal
    beam_area = np.pi*beam_radius**2

    diagonal = (lambda_c + 2.0*mu_c
                + 2.0*shell_thickness/edge_length*shell_normal
                + beam_area/edge_length**2*E_b)
    off_diagonal = (lambda_c
                    + shell_thickness/edge_length*shell_coupling)

    stiffness = np.full((3, 3), off_diagonal)
    np.fill_diagonal(stiffness, diagonal)
    return stiffness


def compute_residuals(reduced_stiffness, analytical_stiffness, hill_mandel,
                      weighted_stiffness):
    """Return the names and values of all verification residuals."""
    if reduced_stiffness.shape != (6, 6):
        raise ValueError("Cubic symmetry requires a 6x6 stiffness matrix.")

    cubic = np.zeros((6, 6))
    normal_diagonal = np.mean(np.diag(reduced_stiffness)[:3])
    normal_coupling = np.mean(
        reduced_stiffness[:3, :3][~np.eye(3, dtype=bool)])
    shear_diagonal = np.mean(np.diag(reduced_stiffness)[3:])

    cubic[:3, :3] = normal_coupling
    np.fill_diagonal(cubic[:3, :3], normal_diagonal)
    np.fill_diagonal(cubic[3:, 3:], shear_diagonal)

    stiffness_norm = np.linalg.norm(reduced_stiffness)
    normal_error = (
        np.linalg.norm(reduced_stiffness[:3, :3] - analytical_stiffness)
        / np.linalg.norm(analytical_stiffness))
    major_symmetry = (
        np.max(np.abs(reduced_stiffness - reduced_stiffness.T))
        / np.max(np.abs(reduced_stiffness)))
    cubic_symmetry = (
        np.linalg.norm(reduced_stiffness - cubic)/stiffness_norm)
    reduction = (
        np.linalg.norm(reduced_stiffness - weighted_stiffness)
        / stiffness_norm)

    names = (
        "normal_block_error",
        "major_symmetry_residual",
        "cubic_symmetry_residual",
        "reduction_residual",
        "max_hill_mandel_residual",
    )
    values = np.array([
        normal_error,
        major_symmetry,
        cubic_symmetry,
        reduction,
        np.max(hill_mandel),
    ])
    return names, values


def assign_constituent_sections(cell, properties, shell_thickness,
                                beam_radius):
    """Create the constituent materials and assign their sections."""
    for name, (youngs_modulus, poissons_ratio) in properties.items():
        material = cell.model.Material(name=name)
        material.Elastic(table=[[youngs_modulus, poissons_ratio]])

    cell.assign_solid_section("CONTINUUM", "CONTINUUM")
    for name in ("SHELL_X", "SHELL_Y", "SHELL_Z"):
        cell.assign_shell_section(name, "SHELL", shell_thickness,
                                  pre_integrate=True)
    for name in ("BEAM_X", "BEAM_Y", "BEAM_Z"):
        cell.assign_beam_section(name, "BEAM", beam_radius,
                                 pre_integrate=True)


def run_stiffness_cases(cell, step_name, strain_magnitude, component_names,
                        collect_diagnostics=False, power_tolerance=None,
                        vtk_directory=None, vtk_element_sets=None):
    """Run six strain cases and return stiffness and optional diagnostics."""
    if len(component_names) != 6:
        raise ValueError("A three-dimensional stiffness requires six cases.")
    if collect_diagnostics and power_tolerance is None:
        raise ValueError(
            "A power tolerance is required for Hill--Mandel diagnostics.")

    stiffness = np.zeros((6, 6))
    hill_mandel = np.zeros(6) if collect_diagnostics else None

    for component, component_name in enumerate(component_names):
        print("\n------------------------------")
        print("Run strain case {} ({}).".format(
            component + 1, component_name))

        def strain_history(time, index=component):
            """Return one ramped macroscopic strain component."""
            strain = np.zeros(6)
            strain[index] = time*strain_magnitude
            return voigt_to_tensor(strain, "STRAIN")

        history = LoadHistory(np.array([0.0, 1.0]), strain_history)
        cell.apply_periodic_bc(step_name, history)
        print("  Run Abaqus analysis.")
        sys.stdout.flush()
        cell.run_analysis()

        print("  Post-process effective stress.")
        sys.stdout.flush()
        cell.open_odb_file()
        try:
            stress = cell.average_stress(step_name)
            stiffness[:, component] = (
                tensor_to_voigt(stress, "STRESS")/strain_magnitude)

            if collect_diagnostics:
                current_time = cell.get_frame_time(step_name, -1)
                previous_time = cell.get_frame_time(step_name, -2)
                strain_increment = (strain_history(current_time)
                                    - strain_history(previous_time))
                macroscopic_work = float(np.sum(stress*strain_increment))
                microscopic_work = cell.average_power_increment(step_name)
                hill_mandel[component] = hill_mandel_residual(
                    macroscopic_work, microscopic_work, power_tolerance)
        finally:
            cell.close_odb_file()

        if vtk_directory is not None:
            print("  Export VTK results.")
            sys.stdout.flush()
            cell.export_to_vtk(
                os.path.join(vtk_directory, component_name),
                field_outputs=["U", "UR"],
                element_sets=vtk_element_sets, clear=True)

    return stiffness, hill_mandel


def export_results(result_directory, reduced_stiffness,
                   analytical_stiffness, hill_mandel,
                   weighted_stiffness):
    """Export the stiffness matrices and verification residuals."""
    if not os.path.isdir(result_directory):
        os.makedirs(result_directory)

    export_to_csv_file(
        reduced_stiffness,
        os.path.join(result_directory, "effective_stiffness_computed"))
    export_to_csv_file(
        reduced_stiffness[:3, :3],
        os.path.join(result_directory, "normal_stiffness_computed"))
    export_to_csv_file(
        analytical_stiffness,
        os.path.join(result_directory, "normal_stiffness_analytical"))
    export_to_csv_file(
        hill_mandel,
        os.path.join(result_directory, "hill_mandel_residuals"),
        header="hill_mandel_residual")
    export_to_csv_file(
        weighted_stiffness,
        os.path.join(result_directory, "effective_stiffness_weighted"))

    residual_names, residuals = compute_residuals(
        reduced_stiffness, analytical_stiffness, hill_mandel,
        weighted_stiffness)
    export_to_csv_file(
        residuals[np.newaxis, :],
        os.path.join(result_directory, "verification_residuals.csv"),
        header=",".join(residual_names))


def main():
    """Run the reduced and ownership-weighted verification."""
    original_directory = os.getcwd()
    try:
        # Material and geometry (all quantities are nondimensional)
        edge_length = 1.0
        shell_thickness = 0.02
        beam_radius = 0.04

        # Values are (Young's modulus, Poisson's ratio)
        properties = {
            "CONTINUUM": (100.0, 0.3),
            "SHELL": (500.0, 0.3),
            "BEAM": (2000.0, 0.3),
        }

        # Mesh, loading, and output settings
        elements_per_edge = 25
        nodes_per_edge = elements_per_edge + 1
        strain_magnitude = 0.01
        power_tolerance = 1.0e-8
        step_name = "LOAD"
        component_names = ("xx", "yy", "zz", "yz", "xz", "xy")

        # Gmsh line, quadrilateral, and hexahedral element types
        element_types = {1: "B31", 3: "S4", 5: "C3D8"}
        vtk_element_sets = [
            "CONTINUUM",
            "SHELL_X", "SHELL_Y", "SHELL_Z",
            "BEAM_X", "BEAM_Y", "BEAM_Z",
        ]

        reduced_directory = os.path.join(here, "abq")
        weighted_directory = os.path.join(here, "abq_weighted")
        result_directory = os.path.join(here, "results")
        vtk_directory = os.path.join(here, "vtk_output")

        # ---------------------------------------------------------------
        # Build the conforming mixed-dimensional mesh once
        model_name = "RVE"
        create_geometry(model_name, edge_length, nodes_per_edge)
        gmsh_cell = GmshRVE(model_name)
        gmsh_cell.add_periodic_mesh_constraints()
        gmsh_cell.generate_mesh()

        # ---------------------------------------------------------------
        # Homogenize the ownership-consistent reduced cell
        if not os.path.isdir(reduced_directory):
            os.makedirs(reduced_directory)
        os.chdir(reduced_directory)

        reduced_cell = AbaqusRVE(gmsh_cell, custom_types=element_types)
        assign_constituent_sections(
            reduced_cell, properties, shell_thickness, beam_radius)
        reduced_cell.model.StaticStep(name=step_name, previous="Initial")

        print("\nHomogenize the reduced cell.")
        reduced_stiffness, hill_mandel = run_stiffness_cases(
            reduced_cell, step_name, strain_magnitude, component_names,
            collect_diagnostics=True, power_tolerance=power_tolerance,
            vtk_directory=os.path.join(vtk_directory, "reduced"),
            vtk_element_sets=vtk_element_sets)

        # ---------------------------------------------------------------
        # Repeat the analysis on the complete ownership-weighted cell
        if not os.path.isdir(weighted_directory):
            os.makedirs(weighted_directory)
        os.chdir(weighted_directory)

        weighted_properties = {
            "CONTINUUM": properties["CONTINUUM"],
            "SHELL": (0.5*properties["SHELL"][0],
                      properties["SHELL"][1]),
            "BEAM": (0.25*properties["BEAM"][0],
                     properties["BEAM"][1]),
        }
        weighted_cell = AbaqusRVE(
            gmsh_cell, custom_types=element_types, remove_secondary=False)
        assign_constituent_sections(
            weighted_cell, weighted_properties, shell_thickness, beam_radius)
        weighted_cell.model.StaticStep(name=step_name, previous="Initial")

        print("\nHomogenize the ownership-weighted complete cell.")
        weighted_stiffness, _ = run_stiffness_cases(
            weighted_cell, step_name, strain_magnitude, component_names,
            vtk_directory=os.path.join(vtk_directory, "weighted"),
            vtk_element_sets=vtk_element_sets + ["SECONDARY_ELEMENTS"])

        # ---------------------------------------------------------------
        # Evaluate and export all verification measures
        analytical_stiffness = analytical_normal_stiffness(
            edge_length, properties["CONTINUUM"], properties["SHELL"],
            properties["BEAM"], shell_thickness, beam_radius)
        print("\nExport verification results.")
        export_results(
            result_directory, reduced_stiffness, analytical_stiffness,
            hill_mandel, weighted_stiffness)
    finally:
        os.chdir(original_directory)
        if gmsh.is_initialized():
            gmsh.finalize()


if __name__ == "__main__":
    # Abaqus replaces stdout; restore the shell that launched the example.
    if sys.__stdout__ is not None:
        sys.stdout = sys.__stdout__
    main()
