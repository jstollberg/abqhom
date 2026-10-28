"""Cyclic small-strain plasticity of a matrix-filled sandwich cell.

The cell combines a soft elastic continuum matrix with perfectly plastic
shell faces and an independently meshed, embedded beam core. The complete
periodic representation contains both shell faces normal to ``X_2``. The
ownership-consistent Abaqus model removes the secondary face before solving
the boundary-value problem.

Run the complete analysis with

    abq2024le cae noGUI=main.py

A regular Python interpreter creates and checks only the Gmsh mesh.
"""

import inspect
import os
import sys

import numpy as np


here = os.path.dirname(os.path.abspath(
    inspect.getframeinfo(inspect.currentframe()).filename))
abqhom_path = os.path.abspath(os.path.join(here, os.pardir,
                                           os.pardir, os.pardir, "src"))
gmsh_path = "C:/Users/jonat/Documents/gmsh/lib"
for path in (abqhom_path, gmsh_path):
    if path not in sys.path:
        sys.path.insert(0, path)

import gmsh

from abqhom.gmshmodel import GmshRVE
from sandwich import create_geometry, nodes_in_physical_group


def cyclic_strain_path(time, reversal_times, reversal_strains):
    """Return the prescribed, fully reversed in-plane strain history."""
    strain = np.zeros((3, 3))
    strain[0, 0] = np.interp(time, reversal_times, reversal_strains)
    return strain


def assign_materials_and_sections(unit_cell, parameters):
    """Assign the elastic matrix and perfectly plastic structural phases."""
    matrix = unit_cell.model.Material(name="MATRIX")
    matrix.Elastic(table=((parameters["matrix_E"],
                           parameters["matrix_nu"]),))

    shell = unit_cell.model.Material(name="SHELL")
    shell.Elastic(table=((parameters["shell_E"],
                          parameters["shell_nu"]),))
    shell.Plastic(table=((parameters["shell_yield"], 0.0),))

    beam = unit_cell.model.Material(name="BEAM")
    beam.Elastic(table=((parameters["beam_E"],
                         parameters["beam_nu"]),))
    beam.Plastic(table=((parameters["beam_yield"], 0.0),))

    unit_cell.assign_solid_section("CONTINUUM_MATRIX", "MATRIX")
    unit_cell.assign_shell_section(
        "SHELL_FACES", "SHELL", parameters["shell_thickness"],
        num_integration_points=5, pre_integrate=False)
    unit_cell.assign_beam_section(
        "BEAM_STRUTS", "BEAM", parameters["beam_radius"],
        pre_integrate=False)


def evaluate_response(unit_cell, step_name, load_path,
                      power_tolerance=1.0e-8):
    """Evaluate effective stress and Hill--Mandel residual per frame."""
    unit_cell.open_odb_file()
    try:
        frame_times = [frame.frameValue for frame in
                       unit_cell.odb.steps[step_name].frames]
        records = []
        previous_strain = None
        for frame, time in enumerate(frame_times):
            strain = load_path(time)
            stress = unit_cell.average_stress(step_name, frame=frame)
            scale = max(np.linalg.norm(stress), 1.0e-30)
            asymmetry = np.linalg.norm(stress - stress.T)/scale

            residual = np.nan
            if previous_strain is not None:
                increment = strain - previous_strain
                macroscopic_work = float(np.sum(stress*increment))
                microscopic_work = unit_cell.average_virtual_power(
                    step_name, frame=frame)
                denominator = max(abs(macroscopic_work),
                                  abs(microscopic_work), power_tolerance)
                residual = abs(macroscopic_work - microscopic_work)/denominator

            records.append(np.concatenate((
                [time, strain[0, 0]], stress.reshape(-1),
                [asymmetry, residual])))
            previous_strain = strain
        return np.asarray(records)
    finally:
        unit_cell.close_odb_file()


def evaluate_phase_plasticity(unit_cell, step_name, phase_sets):
    """Return the maximum equivalent plastic strain in each phase."""
    unit_cell.open_odb_file()
    try:
        step = unit_cell.odb.steps[step_name]
        instance = unit_cell.odb.rootAssembly.instances[unit_cell.model_name]
        regions = [instance.elementSets["ELSET_" + name]
                   for name in phase_sets]

        records = []
        for frame in step.frames:
            peeq = frame.fieldOutputs["PEEQ"]
            maxima = []
            for region in regions:
                values = peeq.getSubset(region=region).values
                maxima.append(max(float(value.data) for value in values))
            records.append([frame.frameValue] + maxima)
    finally:
        unit_cell.close_odb_file()

    return np.asarray(records)


def write_csv(path, data, header):
    """Write numerical results without requiring a plotting dependency."""
    np.savetxt(path, data, delimiter=",", header=header, comments="")


# ---------------------------------------------------------------------------
# Geometry and mesh. The coarse mesh stays below the 1000-node limit of
# Abaqus 2024 Learning Edition, including a margin for reference points.
model_name = "EXAMPLE_3"
edge_length = 1.0
matrix_element_size = edge_length/8.0
elements_per_strut = 16
learning_edition_node_limit = 1000

create_geometry(model_name, edge_length, elements_per_strut)

gmsh_unit_cell = GmshRVE(model_name)
gmsh_unit_cell.add_periodic_mesh_constraints()
gmsh_unit_cell.generate_mesh(
    element_order=1, algorithm_2d=6, algorithm_3d=1,
    element_size=matrix_element_size)
gmsh_unit_cell.build_node_sets()
gmsh_unit_cell.build_element_sets()

matrix_nodes = nodes_in_physical_group(3, "CONTINUUM_MATRIX")
shell_nodes = nodes_in_physical_group(2, "SHELL_FACES")
beam_nodes = nodes_in_physical_group(1, "BEAM_STRUTS")
all_nodes = set(int(node)
                for node in gmsh_unit_cell.node_sets["ALL_NODES"])

if not shell_nodes <= matrix_nodes:
    raise RuntimeError(
        "The shell faces are not conforming with the continuum matrix.")
shared_beam_nodes = beam_nodes & matrix_nodes
if shared_beam_nodes != beam_nodes & shell_nodes:
    raise RuntimeError(
        "Beam nodes shared with the matrix are not the shell corner nodes.")
if len(shared_beam_nodes) != 8:
    raise RuntimeError(
        "The beam core must share exactly its eight corner nodes with the "
        "continuum and shell meshes.")
if len(all_nodes) >= learning_edition_node_limit:
    raise RuntimeError(
        "The mesh contains {} nodes; the Abaqus Learning Edition model must "
        "remain below {}.".format(len(all_nodes),
                                  learning_edition_node_limit))

mesh_file = os.path.join(here, "example_3.msh")
gmsh.write(mesh_file)

print(gmsh_unit_cell)
print("Continuum nodes: {}".format(len(matrix_nodes)))
print("Shell nodes: {}".format(len(shell_nodes)))
print("Beam nodes: {} ({} shared corner nodes)".format(
    len(beam_nodes), len(shared_beam_nodes)))
print("Total nodes: {} / {}".format(
    len(all_nodes), learning_edition_node_limit))
for name in ("CONTINUUM_MATRIX", "SHELL_FACES", "BEAM_STRUTS",
             "SECONDARY_ELEMENTS"):
    print("{} elements: {}".format(
        name, len(gmsh_unit_cell.element_sets[name])))
print("Periodic Gmsh mesh written to {}.".format(mesh_file))


# Stop after mesh generation outside the Abaqus Python environment.
try:
    from abaqusConstants import OFF
except ImportError:
    gmsh.finalize()
    print("Abaqus modules are unavailable; mesh generation is complete.")
    sys.exit(0)

from abqhom.abqmodel import AbaqusRVE, LoadHistory


# ---------------------------------------------------------------------------
# Constituent parameters. Shells and beams have the same elastic constants
# but distinct yield stresses. A single plastic-data point defines ideal J2
# plasticity without hardening.
parameters = {
    "matrix_E": 1.0,
    "matrix_nu": 0.3,
    "shell_E": 100.0,
    "shell_nu": 0.3,
    "shell_yield": 0.4,
    "beam_E": 100.0,
    "beam_nu": 0.3,
    "beam_yield": 0.6,
    "shell_thickness": 0.02,
    "beam_radius": 0.05,
}

# Two fully reversed cycles with a constant strain-rate magnitude. All
# transverse macroscopic strain components remain zero.
reversal_strains = np.asarray([
    0.0, 0.03, -0.03, 0.03, -0.03, 0.03, 0.0])
reversal_times = np.asarray([
    0.0, 0.10, 0.30, 0.50, 0.70, 0.90, 1.00])


def load_path(time):
    return cyclic_strain_path(time, reversal_times, reversal_strains)


history = LoadHistory(reversal_times, load_path)

analysis_directory = os.path.join(here, "abq")
if not os.path.isdir(analysis_directory):
    os.makedirs(analysis_directory)
os.chdir(analysis_directory)

unit_cell = AbaqusRVE(
    gmsh_unit_cell, custom_types={1: "B31"}, remove_secondary=True)
assign_materials_and_sections(unit_cell, parameters)

# The beam endpoints are shared with the matrix and shell meshes. The
# complete beam core can therefore be used as the embedded-element region;
# the shared corner translations already coincide with the host nodes.
beam_region = unit_cell.instance.sets["ELSET_BEAM_STRUTS"]
matrix_region = unit_cell.instance.sets["ELSET_CONTINUUM_MATRIX"]

step_name = "CYCLIC_LOADING"
unit_cell.model.StaticStep(
    name=step_name, previous="Initial", nlgeom=OFF,
    timePeriod=float(reversal_times[-1]), maxNumInc=500,
    initialInc=0.002, minInc=0.002, maxInc=0.002)
# Establish periodicity on the shared boundary nodes before registering the
# embedded constraint for the complete beam core.
unit_cell.apply_periodic_BC(step_name, history)
unit_cell.model.EmbeddedRegion(
    "EMBEDDED_BEAM_CORE", beam_region, matrix_region)
unit_cell.run_analysis(
    field_outputs=["S", "PE", "PEEQ"],
    history_outputs=["ALLIE", "ALLPD", "ALLWK"],
    numCpus=1)

response = evaluate_response(unit_cell, step_name, load_path)
phase_sets = ["SHELL_FACES", "BEAM_STRUTS"]
phase_plasticity = evaluate_phase_plasticity(
    unit_cell, step_name, phase_sets)
result_directory = os.path.join(here, "results")
if not os.path.isdir(result_directory):
    os.makedirs(result_directory)

header = ("time,epsilon_11,sigma_11,sigma_12,sigma_13,"
          "sigma_21,sigma_22,sigma_23,sigma_31,sigma_32,sigma_33,"
          "stress_asymmetry,hill_mandel_residual")
write_csv(os.path.join(result_directory, "cyclic_response.csv"),
          response, header)
write_csv(os.path.join(result_directory, "phase_plasticity.csv"),
          phase_plasticity,
          "time,max_peeq_shell,max_peeq_beam")

frame_times = response[:, 0]
print("\nCyclic effective response:")
print("  reversal | pseudo-time | epsilon_11 | sigma_11")
for reversal, (strain, time) in enumerate(
        zip(reversal_strains[1:-1], reversal_times[1:-1]), start=1):
    frame = int(np.argmin(np.abs(frame_times - time)))
    print("  {:8d} | {:11.2f} | {:+.4f}    | {:+.6e}".format(
        reversal, time, strain, response[frame, 2]))
print("  maximum Hill--Mandel residual = {:.3e}".format(
    np.nanmax(response[:, -1])))
print("  maximum shell PEEQ = {:.6e}".format(
    np.max(phase_plasticity[:, 1])))
print("  maximum beam PEEQ = {:.6e}".format(
    np.max(phase_plasticity[:, 2])))

gmsh.finalize()
print("*** Finished ***")
