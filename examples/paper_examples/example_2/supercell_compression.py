"""Finite-strain compression of a 2 x 2 x 2 Octet supercell.

Two imperfections are compared: the lowest (primitive-periodic) supercell
mode and the lowest genuinely non-primitive supercell mode.
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
paraqus_path = "C:/Users/jonat/Documents/Paraqus/src"
for path in (abqhom_path, gmsh_path, paraqus_path):
    if path not in sys.path:
        sys.path.insert(0, path)

import gmsh

from abqhom.abqmodel import AbaqusRVE, LoadHistory
from abqhom.gmshmodel import GmshRVE
from abqhom.utils import export_to_csv_file
from supercell import create_geometry, apply_eigenmode_imperfection

from abaqusConstants import (ISOTROPIC, NEO_HOOKE, OFF, ON, QUASI_STATIC,
                             VOLUMETRIC_DATA)


def assign_section(supercell, radius, c10, d1, density):
    """Assign a compressible Neo-Hookean circular beam section."""
    material = supercell.model.Material(name="OCTET_MATERIAL")
    material.Hyperelastic(
        materialType=ISOTROPIC, testData=OFF, type=NEO_HOOKE,
        volumetricResponse=VOLUMETRIC_DATA, table=((c10, d1),))
    material.Density(table=((density,),))
    supercell.assign_beam_section(
        "OCTET_BEAMS", "OCTET_MATERIAL", radius, pre_integrate=False)


def frame_times(supercell, step_name):
    """Return all output-frame times."""
    supercell.open_odb_file()
    try:
        return [frame.frameValue
                for frame in supercell.odb.steps[step_name].frames]
    finally:
        supercell.close_odb_file()


def energy_ratios(supercell, step_name):
    """Return peak, final, and global kinetic-energy ratios."""
    supercell.open_odb_file()
    try:
        step = supercell.odb.steps[step_name]
        for region in step.historyRegions.values():
            outputs = region.historyOutputs
            if "ALLIE" not in outputs or "ALLKE" not in outputs:
                continue
            internal = dict(outputs["ALLIE"].data)
            kinetic = dict(outputs["ALLKE"].data)
            records = [(time, internal[time], kinetic[time])
                       for time in sorted(set(internal) & set(kinetic))]
            threshold = 0.02*records[-1][0]
            ratios = [(time, abs(ke)/max(abs(ie), 1.0e-30))
                      for time, ie, ke in records if time >= threshold]
            peak_time, peak_ratio = max(ratios,
                                        key=lambda item: item[1])
            final_ratio = ratios[-1][1]
            global_ratio = (max(abs(ke) for _, _, ke in records)
                            / max(abs(ie) for _, ie, _ in records))
            return peak_ratio, peak_time, final_ratio, global_ratio
    finally:
        supercell.close_odb_file()
    raise RuntimeError("ALLIE and ALLKE are unavailable in the ODB.")


def evaluate_response(supercell, step_name, load_path):
    """Evaluate the nominal and Cauchy axial stress in all frames."""
    response = []
    supercell.open_odb_file()
    try:
        step = supercell.odb.steps[step_name]
        instance = supercell.odb.rootAssembly.instances[
            supercell.model_name]
        node_coordinates = {
            node.label: np.asarray(node.coordinates, dtype=float)
            for node in instance.nodes}
        volume = supercell.get_volume()

        for frame in step.frames:
            gradient = load_path(frame.frameValue)
            deformation_gradient = np.eye(3) + gradient
            outputs = frame.fieldOutputs
            forces = {}
            for component in range(3):
                for value in outputs[
                        "NFORC{}".format(component + 1)].values:
                    if value.elementLabel is None:
                        continue
                    key = (value.elementLabel, value.nodeLabel)
                    if key not in forces:
                        forces[key] = np.zeros(3)
                    forces[key][component] = value.data

            nominal_stress = np.zeros((3, 3))
            for (_, node_label), force in forces.items():
                nominal_stress += np.outer(
                    force, node_coordinates[node_label])
            nominal_stress *= -1.0/volume
            jacobian = np.linalg.det(deformation_gradient)
            cauchy_stress = nominal_stress.dot(
                deformation_gradient.T)/jacobian
            response.append([
                frame.frameValue, 1.0 + gradient[0, 0],
                nominal_stress[0, 0], cauchy_stress[0, 0]])
    finally:
        supercell.close_odb_file()
    return np.asarray(response)


def run_case(gmsh_supercell, case_name, mode_number, parameters):
    """Solve one compression case initialized by a selected eigenmode."""
    mode_directory = os.path.join(
        here, "supercell_results",
        "buckling_2x2x2_n{}".format(parameters["elements_per_strut"]))
    residual = apply_eigenmode_imperfection(
        gmsh_supercell, mode_directory, mode_number,
        parameters["imperfection_amplitude"])

    analysis_directory = os.path.join(here, "sc_abq", case_name)
    if not os.path.isdir(analysis_directory):
        os.makedirs(analysis_directory)
    os.chdir(analysis_directory)

    supercell = AbaqusRVE(
        gmsh_supercell, custom_types={1: "B31"}, remove_secondary=True)
    assign_section(
        supercell, parameters["radius"], parameters["c10"],
        parameters["d1"], parameters["density"])

    duration = parameters["step_duration"]
    final_gradient = parameters["final_gradient"]

    def load_path(time):
        gradient = np.zeros((3, 3))
        gradient[0, 0] = final_gradient*time/duration
        return gradient

    step_name = "LOAD"
    supercell.model.ImplicitDynamicsStep(
        name=step_name, previous="Initial", application=QUASI_STATIC,
        nlgeom=ON, timePeriod=duration, maxNumInc=8000,
        initialInc=parameters["time_increment"],
        minInc=1.0e-7*duration, maxInc=parameters["time_increment"])
    history = LoadHistory(np.asarray([0.0, duration]), load_path)
    supercell.apply_periodic_BC(step_name, history)
    supercell.run_analysis(
        field_outputs=["S"],
        history_outputs=["ALLIE", "ALLKE", "ALLVD", "ALLWK"])

    times = frame_times(supercell, step_name)
    if not np.isclose(times[-1], duration):
        raise RuntimeError(
            "The {} case stopped at time {}; expected {}.".format(
                case_name, times[-1], duration))

    response = evaluate_response(supercell, step_name, load_path)
    energies = energy_ratios(supercell, step_name)
    result_directory = os.path.join(
        here, "supercell_results", case_name)
    if not os.path.isdir(result_directory):
        os.makedirs(result_directory)
    export_to_csv_file(
        response, os.path.join(result_directory, "stress_curve"))
    export_to_csv_file(
        np.asarray([[mode_number, residual] + list(energies)]),
        os.path.join(result_directory, "verification_metrics"))

    print("\n{}:".format(case_name))
    print("  eigenmode              : {}".format(mode_number))
    print("  final stretch          : {:.8e}".format(response[-1, 1]))
    print("  final Cauchy stress    : {:.8e}".format(response[-1, 3]))
    print("  peak ALLKE/ALLIE       : {:.3e}".format(energies[0]))
    return response


model_name = "SC2_COMP"
edge_length = 1.0
cells_per_axis = 2
elements_per_strut = int(os.environ.get("OCTET_ELEMENTS_PER_STRUT", "16"))
beam_radius = 0.04
youngs_modulus = 100.0
poissons_ratio = 0.3
step_duration = float(os.environ.get("OCTET_STEP_DURATION", "1.0"))

create_geometry(
    model_name, edge_length, elements_per_strut, cells_per_axis)
gmsh_supercell = GmshRVE(model_name)
gmsh_supercell.add_periodic_mesh_constraints()
gmsh_supercell.generate_mesh()
gmsh_supercell.build_node_sets()
gmsh_supercell.build_element_sets()

original_nodes = {}
for node_tag in gmsh_supercell.node_sets["ALL_NODES"]:
    coordinates, parametric, _, _ = gmsh.model.mesh.get_node(int(node_tag))
    original_nodes[int(node_tag)] = (
        np.asarray(coordinates, dtype=float), parametric)

parameters = {
    "elements_per_strut": elements_per_strut,
    "radius": beam_radius,
    "imperfection_amplitude": 0.01*beam_radius,
    "c10": youngs_modulus/(4.0*(1.0 + poissons_ratio)),
    "d1": 6.0*(1.0 - 2.0*poissons_ratio)/youngs_modulus,
    "density": 2.5e-7,
    "step_duration": step_duration,
    "time_increment": 0.0025*step_duration,
    "final_gradient": -0.40,
}

case_tag = "n{}".format(elements_per_strut)
cases = (("comp_{}_mode_01".format(case_tag), 1),
         ("comp_{}_mode_03".format(case_tag), 3))
selected_mode = os.environ.get("OCTET_SUPERCELL_MODE")
for case_name, mode_number in cases:
    if (selected_mode is not None
            and mode_number != int(selected_mode)):
        continue
    for node_tag, (coordinates, parametric) in original_nodes.items():
        gmsh.model.mesh.set_node(node_tag, coordinates, parametric)
    run_case(gmsh_supercell, case_name, mode_number, parameters)

os.chdir(here)
gmsh.finalize()
print("*** Finished 2 x 2 x 2 Octet compression comparison ***")
