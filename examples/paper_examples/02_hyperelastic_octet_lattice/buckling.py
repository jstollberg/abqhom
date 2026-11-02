# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-or-later

"""Compute periodic buckling modes of the Octet cell and supercell.

Run this script before ``main.py``. It writes the normalized translational
eigenmodes used to construct the geometric imperfections of both cell sizes.
"""

import os
import re

import numpy as np

here = os.getcwd()

import gmsh
from abaqusConstants import (BUCKLING_MODES, STRESS_PERTURBATION,
                             SUBSPACE, UNSET)

from abqhom.abqmodel import AbaqusRVE, LoadHistory
from abqhom.constants import PERIODIC_PAIRS
from abqhom.gmshmodel import GmshRVE
from abqhom.utils import export_to_csv_file
from octet import create_geometry


def remove_temporary_loading(cell, setup_step):
    """Retain the periodic equations and reference points only."""
    model = cell.model
    for name in list(model.boundaryConditions.keys()):
        if name.startswith("DIRICHLET"):
            del model.boundaryConditions[name]
    for name in list(model.amplitudes.keys()):
        if name.startswith("DIRICHLET"):
            del model.amplitudes[name]
    del model.steps[setup_step]


def add_buckling_boundary_conditions(cell, step_name,
                                     compression_gradient):
    """Define the perturbation load and periodic buckling-mode space."""
    model = cell.model
    assembly = model.rootAssembly
    cell_sets = cell.instance.sets

    # The stress-perturbation anchor follows the affine compression, whereas
    # the eigenmode anchor suppresses only the arbitrary rigid translation.
    anchor_region, anchor = cell.select_rigid_body_anchor()
    anchor_displacement = compression_gradient.dot(
        np.asarray(anchor.coordinates, dtype=float))
    model.DisplacementBC(
        name="BUCKLE_STRESS_ANCHOR", createStepName=step_name,
        region=anchor_region, amplitude=UNSET,
        u1=anchor_displacement[0], u2=anchor_displacement[1],
        u3=anchor_displacement[2], buckleCase=STRESS_PERTURBATION)
    model.DisplacementBC(
        name="BUCKLE_MODE_ANCHOR", createStepName=step_name,
        region=anchor_region, amplitude=UNSET,
        u1=0.0, u2=0.0, u3=0.0, buckleCase=BUCKLING_MODES)

    for secondary, primary in PERIODIC_PAIRS.items():
        secondary_set = "NSET_" + secondary
        primary_set = "NSET_" + primary
        if primary_set not in cell_sets or secondary_set not in cell_sets:
            continue

        point_name = "REFERENCE_POINT_{}_{}".format(primary, secondary)
        region = assembly.instances[point_name].sets["REFERENCE_POINT"]
        x_primary = np.asarray(
            cell_sets[primary_set].nodes[0].coordinates, dtype=float)
        x_secondary = np.asarray(
            cell_sets[secondary_set].nodes[0].coordinates, dtype=float)
        displacement = compression_gradient.dot(x_secondary - x_primary)

        # Reference-point translations generate the pre-buckling stress
        # pattern but remain fixed in the homogeneous eigenmode problem.
        model.DisplacementBC(
            name="BUCKLE_STRESS_{}".format(point_name),
            createStepName=step_name, region=region, amplitude=UNSET,
            u1=displacement[0], u2=displacement[1],
            u3=displacement[2], buckleCase=STRESS_PERTURBATION)
        model.DisplacementBC(
            name="BUCKLE_MODE_{}".format(point_name),
            createStepName=step_name, region=region, amplitude=UNSET,
            u1=0.0, u2=0.0, u3=0.0, buckleCase=BUCKLING_MODES)


def periodicity_residual(cell, displacement):
    """Return the largest normalized primary--secondary mode mismatch."""
    maximum = max(np.linalg.norm(value)
                  for value in displacement.values())
    if maximum <= 0.0:
        raise RuntimeError("The buckling mode has zero displacement.")

    residual = 0.0
    sets = cell.instance.sets
    for secondary, primary in PERIODIC_PAIRS.items():
        secondary_set = "NSET_" + secondary
        primary_set = "NSET_" + primary
        pair_available = primary_set in sets and secondary_set in sets
        if not pair_available:
            continue

        primary_nodes = sets[primary_set].nodes
        secondary_nodes = sets[secondary_set].nodes
        if len(primary_nodes) != len(secondary_nodes):
            raise RuntimeError(
                "Periodic node sets '{}' and '{}' have different sizes."
                .format(primary_set, secondary_set))

        for primary_node, secondary_node in zip(
                primary_nodes, secondary_nodes):
            primary_value = displacement[primary_node.label]
            secondary_value = displacement[secondary_node.label]
            mismatch = np.linalg.norm(secondary_value - primary_value)
            residual = max(residual, mismatch)
    return residual / maximum


def primitive_periodicity_residual(instance, displacement, edge_length):
    """Measure whether a supercell mode repeats after one primitive cell."""
    coordinates = {
        node.label: np.asarray(node.coordinates, dtype=float)
        for node in instance.nodes}
    coordinate_to_label = {
        tuple(np.round(value, 8)): label
        for label, value in coordinates.items()}
    maximum = max(np.linalg.norm(value)
                  for value in displacement.values())
    if maximum <= 0.0:
        raise RuntimeError("The buckling mode has zero displacement.")

    residual = 0.0
    comparisons = 0
    for label, coordinate in coordinates.items():
        for axis in range(3):
            shifted = coordinate.copy()
            shifted[axis] += edge_length
            other_label = coordinate_to_label.get(
                tuple(np.round(shifted, 8)))
            if other_label is None:
                continue
            mismatch = np.linalg.norm(
                displacement[other_label] - displacement[label])
            residual = max(residual, mismatch)
            comparisons += 1

    if comparisons == 0:
        raise RuntimeError(
            "No node pairs separated by one primitive cell were found.")
    return residual / maximum, comparisons


def read_eigenvalue(frame):
    """Read the eigenvalue stored in an Abaqus frame description."""
    match = re.search(
        r"EigenValue\s*=\s*([+\-0-9.Ee]+)", frame.description,
        flags=re.IGNORECASE)
    if match is None:
        raise RuntimeError(
            "Could not read the eigenvalue from frame description '{}'."
            .format(frame.description))
    return float(match.group(1))


def export_modes(cell, step_name, result_directory,
                 primitive_edge_length=None):
    """Export normalized modes, eigenvalues, and periodicity residuals."""
    if not os.path.isdir(result_directory):
        os.makedirs(result_directory)

    cell.open_odb_file()
    try:
        step = cell.odb.steps[step_name]
        instance = cell.odb.rootAssembly.instances[cell.model_name]
        records = []
        for mode_number, frame in enumerate(step.frames[1:], start=1):
            values = frame.fieldOutputs["U"].getSubset(
                region=instance).values
            displacement = {
                value.nodeLabel: np.asarray(value.data, dtype=float)
                for value in values}
            maximum = max(np.linalg.norm(value)
                          for value in displacement.values())
            if maximum <= 0.0:
                raise RuntimeError(
                    "Buckling mode {} has zero displacement."
                    .format(mode_number))

            # Reference coordinates make the exported mode self-validating
            # when it is transferred to a separately regenerated Gmsh mesh.
            coordinates = {
                node.label: np.asarray(node.coordinates, dtype=float)
                for node in instance.nodes}
            mode_data = np.asarray([
                [label] + list(coordinates[label]) + list(value / maximum)
                for label, value in sorted(displacement.items())])
            export_to_csv_file(
                mode_data, os.path.join(
                    result_directory, "mode_{:02d}".format(mode_number)))

            eigenvalue = read_eigenvalue(frame)
            residual = periodicity_residual(cell, displacement)
            record = [mode_number, eigenvalue, residual]
            message = ("mode {:2d}: eigenvalue = {:.8e}, outer residual = "
                       "{:.3e}".format(
                           mode_number, eigenvalue, residual))
            if primitive_edge_length is not None:
                primitive_residual, comparisons = (
                    primitive_periodicity_residual(
                        instance, displacement, primitive_edge_length))
                record.extend([primitive_residual, comparisons])
                message += ", primitive residual = {:.3e}".format(
                    primitive_residual)
            records.append(record)
            print(message)
    finally:
        cell.close_odb_file()

    records = np.asarray(records)
    export_to_csv_file(
        records, os.path.join(result_directory, "eigenvalues"))
    return records


def run_buckling(model_name, cells_per_axis, elements_per_strut,
                 number_of_modes, analysis_directory, result_directory,
                 vtk_directory, max_iterations=None):
    """Run one periodic Octet eigenvalue-buckling analysis."""
    original_directory = os.getcwd()
    try:
        # Material and geometry (all quantities are nondimensional)
        edge_length = 1.0
        beam_radius = 0.04
        youngs_modulus = 100.0
        poissons_ratio = 0.3

        setup_step = "PERIODIC_SETUP"
        step_name = "BUCKLING"

        # ---------------------------------------------------------------
        # Build the perfect periodic mesh
        create_geometry(
            model_name, edge_length, elements_per_strut, cells_per_axis)
        gmsh_cell = GmshRVE(model_name)
        gmsh_cell.add_periodic_mesh_constraints()
        gmsh_cell.generate_mesh()

        # ---------------------------------------------------------------
        # Construct the reduced linear model
        if not os.path.isdir(analysis_directory):
            os.makedirs(analysis_directory)
        os.chdir(analysis_directory)

        cell = AbaqusRVE(gmsh_cell, custom_types={1: "B31"})
        material = cell.model.Material(name="OCTET_MATERIAL")
        material.Elastic(table=((youngs_modulus, poissons_ratio),))
        cell.assign_beam_section(
            "OCTET_BEAMS", "OCTET_MATERIAL", beam_radius,
            pre_integrate=True)

        # Generate the periodic equations and auxiliary reference points
        # through the same boundary-condition implementation as the
        # finite-strain analyses. The temporary zero loading is discarded.
        cell.model.StaticStep(name=setup_step, previous="Initial")

        def zero_gradient(_time):
            """Return a zero macroscopic displacement gradient."""
            return np.zeros((3, 3))

        zero_history = LoadHistory(np.asarray([0.0, 1.0]), zero_gradient)
        cell.apply_periodic_BC(setup_step, zero_history)
        remove_temporary_loading(cell, setup_step)

        # ---------------------------------------------------------------
        # Compute the periodic modes under unit macroscopic compression
        step_options = {
            "name": step_name,
            "previous": "Initial",
            "numEigen": number_of_modes,
            "eigensolver": SUBSPACE,
        }
        if max_iterations is not None:
            step_options["maxIterations"] = max_iterations
        cell.model.BuckleStep(**step_options)
        compression_gradient = np.zeros((3, 3))
        compression_gradient[0, 0] = -1.0
        add_buckling_boundary_conditions(
            cell, step_name, compression_gradient)
        cell.run_analysis(scratch=analysis_directory)

        primitive_edge_length = (
            edge_length if cells_per_axis > 1 else None)
        modes = export_modes(
            cell, step_name, result_directory, primitive_edge_length)
        if len(modes) != number_of_modes:
            raise RuntimeError(
                "Expected {} modes, but Abaqus returned {}."
                .format(number_of_modes, len(modes)))

        cell.export_to_vtk(
            vtk_directory, field_outputs=["U", "UR"],
            element_sets=["OCTET_BEAMS"], clear=True)

        print("*** Finished periodic {} x {} x {} Octet buckling analysis "
              "***".format(
                  cells_per_axis, cells_per_axis, cells_per_axis))
    finally:
        os.chdir(original_directory)
        if gmsh.is_initialized():
            gmsh.finalize()


def main(include_supercell=True):
    """Run the primitive-cell and optional supercell analyses."""
    elements_per_strut = 16
    run_buckling(
        model_name="EXAMPLE_2_OCTET_BUCKLING",
        cells_per_axis=1,
        elements_per_strut=elements_per_strut,
        number_of_modes=2,
        analysis_directory=os.path.join(here, "abq", "buckling"),
        result_directory=os.path.join(here, "results", "buckling"),
        vtk_directory=os.path.join(here, "vtk_output", "buckling"))

    if include_supercell:
        result_tag = "buckling_2x2x2_n{}".format(elements_per_strut)
        run_buckling(
            model_name="SC2_BUCKLING",
            cells_per_axis=2,
            elements_per_strut=elements_per_strut,
            number_of_modes=2,
            analysis_directory=os.path.join(
                here, "supercell_abq",
                "n{}".format(elements_per_strut)),
            result_directory=os.path.join(
                here, "supercell_results", result_tag),
            vtk_directory=os.path.join(
                here, "supercell_vtk_output", result_tag),
            # The nearly degenerate lowest eigenspace can require more
            # subspace iterations than the primitive-cell problem.
            max_iterations=100)


if __name__ == "__main__":
    main(include_supercell=True)
