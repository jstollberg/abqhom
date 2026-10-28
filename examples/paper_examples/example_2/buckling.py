"""Periodic eigenvalue-buckling analysis for the perfect Octet cell.

Run this script before ``main.py`` to generate the periodic imperfection
mode used by the finite-strain analyses.
"""

import inspect
import os
import re
import sys

import numpy as np


here = os.path.dirname(os.path.abspath(
    inspect.getframeinfo(inspect.currentframe()).filename))
abqhom_path = os.path.abspath(os.path.join(here, os.pardir,
                                           os.pardir, os.pardir, "src"))
gmsh_path = "C:/Users/jonat/Documents/gmsh/lib"
paraqus_path = "C:/Users/jonat/Documents/Paraqus/src"
for path in (abqhom_path, gmsh_path, paraqus_path):
    if path not in sys.path:
        sys.path.insert(0, path)

import gmsh

from abqhom.abqmodel import AbaqusRVE, LoadHistory
from abqhom.constants import PERIODIC_PAIRS
from abqhom.gmshmodel import GmshRVE
from abqhom.utils import export_to_csv_file
from octet import create_geometry

from abaqusConstants import (BUCKLING_MODES, STRESS_PERTURBATION,
                             SUBSPACE, UNSET)


def remove_temporary_boundary_conditions(unit_cell, setup_step):
    """Retain the periodic equations and reference points only."""
    model = unit_cell.model
    for name in list(model.boundaryConditions.keys()):
        if name.startswith("DIRICHLET"):
            del model.boundaryConditions[name]
    for name in list(model.amplitudes.keys()):
        if name.startswith("DIRICHLET"):
            del model.amplitudes[name]
    del model.steps[setup_step]


def add_buckling_boundary_conditions(unit_cell, step_name,
                                     compression_gradient):
    """Define the perturbation load and fixed-macrostrain mode spaces."""
    model = unit_cell.model
    assembly = model.rootAssembly
    rve_sets = unit_cell.instance.sets

    anchor_region, anchor = unit_cell._select_rigid_body_anchor()
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

    for slave, master in PERIODIC_PAIRS.items():
        slave_set = "NSET_" + slave
        master_set = "NSET_" + master
        if master_set not in rve_sets or slave_set not in rve_sets:
            continue

        point_name = "REFERENCE_POINT_{}_{}".format(master, slave)
        region = assembly.instances[point_name].sets["REFERENCE_POINT"]
        x_minus = np.asarray(rve_sets[master_set].nodes[0].coordinates,
                             dtype=float)
        x_plus = np.asarray(rve_sets[slave_set].nodes[0].coordinates,
                            dtype=float)
        displacement = compression_gradient.dot(x_plus - x_minus)

        model.DisplacementBC(
            name="BUCKLE_STRESS_{}".format(point_name),
            createStepName=step_name, region=region, amplitude=UNSET,
            u1=displacement[0], u2=displacement[1],
            u3=displacement[2], buckleCase=STRESS_PERTURBATION)
        model.DisplacementBC(
            name="BUCKLE_MODE_{}".format(point_name),
            createStepName=step_name, region=region, amplitude=UNSET,
            u1=0.0, u2=0.0, u3=0.0, buckleCase=BUCKLING_MODES)


def periodic_residual(unit_cell, displacement):
    """Return the largest normalized slave-master mode mismatch."""
    maximum = max(np.linalg.norm(value)
                  for value in displacement.values())
    residual = 0.0
    for slave, master in PERIODIC_PAIRS.items():
        slave_set = "NSET_" + slave
        master_set = "NSET_" + master
        if (slave_set not in unit_cell.instance.sets
                or master_set not in unit_cell.instance.sets):
            continue
        slave_nodes = unit_cell.instance.sets[slave_set].nodes
        master_nodes = unit_cell.instance.sets[master_set].nodes
        for slave_node, master_node in zip(slave_nodes, master_nodes):
            mismatch = np.linalg.norm(
                displacement[slave_node.label]
                - displacement[master_node.label])
            residual = max(residual, mismatch)
    return residual/max(maximum, 1.0e-30)


def export_modes(unit_cell, step_name, result_directory):
    """Export eigenvalues, normalized modes, and periodicity diagnostics."""
    if not os.path.isdir(result_directory):
        os.makedirs(result_directory)

    unit_cell.open_odb_file()
    try:
        step = unit_cell.odb.steps[step_name]
        instance = unit_cell.odb.rootAssembly.instances[unit_cell.model_name]
        records = []
        for mode_number in range(1, len(step.frames)):
            frame = step.frames[mode_number]
            values = frame.fieldOutputs["U"].getSubset(
                region=instance).values
            displacement = {
                value.nodeLabel: np.asarray(value.data, dtype=float)
                for value in values}
            maximum = max(np.linalg.norm(value)
                          for value in displacement.values())
            mode_data = np.asarray([
                [label] + list(value/maximum)
                for label, value in sorted(displacement.items())])
            export_to_csv_file(
                mode_data, os.path.join(
                    result_directory, "mode_{:02d}".format(mode_number)))

            match = re.search(
                r"EigenValue\s*=\s*([+\-0-9.Ee]+)", frame.description,
                flags=re.IGNORECASE)
            if match is None:
                raise RuntimeError(
                    "Could not read the eigenvalue from frame description "
                    "'{}'.".format(frame.description))
            eigenvalue = float(match.group(1))
            residual = periodic_residual(unit_cell, displacement)
            records.append([mode_number, eigenvalue, residual])
            print("mode {:2d}: eigenvalue = {:.8e}, periodic residual = "
                  "{:.3e}".format(mode_number, eigenvalue, residual))
    finally:
        unit_cell.close_odb_file()

    records = np.asarray(records)
    export_to_csv_file(
        records, os.path.join(result_directory, "eigenvalues"))
    return records


model_name = "EXAMPLE_2_OCTET_BUCKLING"
edge_length = 1.0
elements_per_strut = int(os.environ.get("OCTET_ELEMENTS_PER_STRUT", "16"))
beam_radius = 0.04
youngs_modulus = 100.0
poissons_ratio = 0.3
number_of_modes = 20

create_geometry(model_name, edge_length, elements_per_strut)
gmsh_unit_cell = GmshRVE(model_name)
gmsh_unit_cell.add_periodic_mesh_constraints()
gmsh_unit_cell.generate_mesh()
gmsh_unit_cell.build_node_sets()
gmsh_unit_cell.build_element_sets()

result_tag = os.environ.get("OCTET_RESULT_TAG", "final")
analysis_directory = os.path.join(
    here, result_tag + "_abq", "buckling")
if not os.path.isdir(analysis_directory):
    os.makedirs(analysis_directory)
os.chdir(analysis_directory)

unit_cell = AbaqusRVE(
    gmsh_unit_cell, custom_types={1: "B31"}, remove_secondary=True)
material = unit_cell.model.Material(name="OCTET_MATERIAL")
material.Elastic(table=((youngs_modulus, poissons_ratio),))
unit_cell.assign_beam_section(
    "OCTET_BEAMS", "OCTET_MATERIAL", beam_radius, pre_integrate=True)

# Construct the node-pair equations and reference points through the standard
# periodic-BC implementation. The temporary zero-valued loading data are then
# replaced by the load and mode boundary conditions of the buckling step.
setup_step = "PERIODIC_SETUP"
unit_cell.model.StaticStep(name=setup_step, previous="Initial")
zero_history = LoadHistory(
    np.asarray([0.0, 1.0]), lambda time: np.zeros((3, 3)))
unit_cell.apply_periodic_BC(setup_step, zero_history)
remove_temporary_boundary_conditions(unit_cell, setup_step)

step_name = "BUCKLING"
unit_cell.model.BuckleStep(
    name=step_name, previous="Initial", numEigen=number_of_modes,
    eigensolver=SUBSPACE)
compression_gradient = np.zeros((3, 3))
compression_gradient[0, 0] = -1.0
add_buckling_boundary_conditions(
    unit_cell, step_name, compression_gradient)

unit_cell.run_analysis(scratch=analysis_directory)

result_directory = os.path.join(
    here, result_tag + "_results", "buckling")
results = export_modes(unit_cell, step_name, result_directory)
if len(results) != number_of_modes:
    raise RuntimeError(
        "Expected {} modes, but Abaqus returned {}.".format(
            number_of_modes, len(results)))

unit_cell.export_to_vtk(
    os.path.join(here, result_tag + "_vtk", "buckling"),
    field_outputs=["U", "UR"], element_sets=["OCTET_BEAMS"], clear=True)

os.chdir(here)
gmsh.finalize()
print("*** Finished periodic Octet buckling analysis ***")
