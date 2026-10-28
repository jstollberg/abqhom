"""
Mixed-dimensional affine patch cell for the first paper example.

The cuboidal unit cell contains a conforming hexahedral continuum mesh,
quadrilateral shell meshes on all six faces, and beam meshes on all twelve
edges.  The three constituent dimensions share the same Gmsh nodes.  In the
complete periodic representation, the continuum elements have unit
ownership, the face-supported shell elements have ownership one half, and
the edge-supported beam elements have ownership one quarter.

The six independent macroscopic strain cases provide the complete effective
stiffness matrix and verify the Hill--Mandel condition. Run with

    abaqus cae noGUI=example_1.py

or with a Python interpreter for which the Gmsh module is available.
"""

import inspect
import os
import sys

import numpy as np


# ---------------------------------------------------------------
# Locate this example and put the shipped abqhom source first on the
# path. Under 'abaqus cae noGUI=example_1.py' __file__ is undefined, so the
# script location is obtained via inspect instead.
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
from abqhom.abqmodel import AbaqusRVE, LoadHistory
from abqhom.utils import (export_to_csv_file, hill_mandel_residual,
                          tensor_to_voigt, voigt_to_tensor)


def _add_physical_group(dimension, entities, name):
    """Add and name one physical group."""
    group = gmsh.model.add_physical_group(dimension, entities)
    gmsh.model.set_physical_name(dimension, group, name)


def create_geometry(model_name, edge_length=1.0, nodes_per_edge=26):
    """Create the conforming continuum--shell--beam unit cell."""
    if not gmsh.is_initialized():
        gmsh.initialize()
    if model_name in gmsh.model.list():
        gmsh.model.set_current(model_name)
        gmsh.model.remove()
    gmsh.model.add(model_name)

    gmsh.option.set_number("Geometry.Tolerance", 1.0e-8*edge_length)

    h = 0.5*edge_length
    g = gmsh.model.geo

    # Corner points. The numbering follows the four corners of the face at
    # y=-h and then the corresponding corners of the face at y=+h.
    p1 = g.add_point(-h, -h, -h)
    p2 = g.add_point(-h, -h,  h)
    p3 = g.add_point( h, -h,  h)
    p4 = g.add_point( h, -h, -h)
    p5 = g.add_point(-h,  h, -h)
    p6 = g.add_point(-h,  h,  h)
    p7 = g.add_point( h,  h,  h)
    p8 = g.add_point( h,  h, -h)

    # Edge curves, grouped below according to their axial direction.
    lines = [g.add_line(p1, p2),
             g.add_line(p2, p3),
             g.add_line(p3, p4),
             g.add_line(p4, p1),
             g.add_line(p5, p6),
             g.add_line(p6, p7),
             g.add_line(p7, p8),
             g.add_line(p8, p5),
             g.add_line(p1, p5),
             g.add_line(p2, p6),
             g.add_line(p3, p7),
             g.add_line(p4, p8)]

    # Boundary surfaces. Their order is y-, y+, x-, z+, x+, z-.
    surface_connectivity = [[lines[0], lines[1], lines[2], lines[3]],
                            [lines[4], lines[5], lines[6], lines[7]],
                            [lines[0], lines[9], lines[4], lines[8]],
                            [lines[1], lines[10], lines[5], lines[9]],
                            [lines[2], lines[11], lines[6], lines[10]],
                            [lines[3], lines[8], lines[7], lines[11]]]
    surfaces = []
    for connectivity in surface_connectivity:
        loop = g.add_curve_loop(connectivity, reorient=True)
        surfaces.append(g.add_plane_surface([loop]))

    surface_loop = g.add_surface_loop(surfaces)
    volume = g.add_volume([surface_loop])

    # A structured conforming mesh gives common nodes to the continuum,
    # shell, and beam elements.
    for line in lines:
        g.mesh.set_transfinite_curve(line, nodes_per_edge)
    for surface in surfaces:
        g.mesh.set_transfinite_surface(surface)
        g.mesh.set_recombine(2, surface)
    g.mesh.set_transfinite_volume(volume)

    g.synchronize()

    # The orientation-specific groups allow the surface and line
    # contributions to the effective response to be checked separately.
    _add_physical_group(3, [volume], "CONTINUUM")
    _add_physical_group(2, [surfaces[2], surfaces[4]], "SHELL_X")
    _add_physical_group(2, [surfaces[0], surfaces[1]], "SHELL_Y")
    _add_physical_group(2, [surfaces[3], surfaces[5]], "SHELL_Z")
    _add_physical_group(1, [lines[1], lines[3], lines[5], lines[7]],
                        "BEAM_X")
    _add_physical_group(1, [lines[8], lines[9], lines[10], lines[11]],
                        "BEAM_Y")
    _add_physical_group(1, [lines[0], lines[2], lines[4], lines[6]],
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


def cubic_symmetry_residual(stiffness):
    """Return the relative Frobenius distance from cubic symmetry."""
    if stiffness.shape != (6, 6):
        raise ValueError("Cubic symmetry requires a 6x6 stiffness matrix.")

    cubic = np.zeros((6, 6))
    normal_diagonal = np.mean(np.diag(stiffness)[:3])
    normal_coupling = np.mean(
        stiffness[:3, :3][~np.eye(3, dtype=bool)])
    shear_diagonal = np.mean(np.diag(stiffness)[3:])

    cubic[:3, :3] = normal_coupling
    np.fill_diagonal(cubic[:3, :3], normal_diagonal)
    np.fill_diagonal(cubic[3:, 3:], shear_diagonal)

    scale = np.linalg.norm(stiffness)
    if scale == 0.0:
        return 0.0
    return np.linalg.norm(stiffness - cubic)/scale


def assign_constituent_sections(unit_cell, properties, shell_thickness,
                                beam_radius):
    """Assign the materials and sections of the mixed-dimensional cell."""
    for name, (E, nu) in properties.items():
        material = unit_cell.model.Material(name=name)
        material.Elastic(table=[[E, nu]])

    unit_cell.assign_solid_section("CONTINUUM", "CONTINUUM")
    for name in ("SHELL_X", "SHELL_Y", "SHELL_Z"):
        unit_cell.assign_shell_section(name, "SHELL", shell_thickness,
                                       pre_integrate=True)
    for name in ("BEAM_X", "BEAM_Y", "BEAM_Z"):
        unit_cell.assign_beam_section(name, "BEAM", beam_radius,
                                      pre_integrate=True)


def _get_internal_energy(unit_cell, step_name, frame=-1):
    """Read the total internal energy from the analysis history."""
    unit_cell.open_odb_file()
    try:
        step = unit_cell.odb.steps[step_name]
        energies = []
        for region in step.historyRegions.values():
            if "ALLIE" in region.historyOutputs:
                energies.append(
                    region.historyOutputs["ALLIE"].data[frame][1])
    finally:
        unit_cell.close_odb_file()

    if len(energies) != 1:
        raise RuntimeError(
            "Expected one global ALLIE history output, found {}."
            .format(len(energies)))
    return float(energies[0])


def _relative_residual(first, second, tolerance=1.0e-14):
    """Return a symmetric normalized difference."""
    scale = max(abs(first), abs(second), tolerance)
    return abs(first - second)/scale


# ---------------------------------------------------------------
# Geometry, discretization, and constituent properties. The unit system is
# arbitrary but must be used consistently.
model_name = "EXAMPLE_1"
edge_length = 1.0
nodes_per_edge = 26
shell_thickness = 0.02
beam_radius = 0.04
properties = {"CONTINUUM": (100.0, 0.3),
              "SHELL":     (500.0, 0.3),
              "BEAM":      (2000.0, 0.3)}

create_geometry(model_name, edge_length, nodes_per_edge)

unit_cell = GmshRVE(model_name)
unit_cell.add_periodic_mesh_constraints()
unit_cell.generate_mesh()
unit_cell.build_node_sets()
unit_cell.build_element_sets()

print("\nMesh summary:")
print("  nodes                 : {}".format(
    len(unit_cell.node_sets["ALL_NODES"])))
for name in ("CONTINUUM", "SHELL_X", "SHELL_Y", "SHELL_Z",
             "BEAM_X", "BEAM_Y", "BEAM_Z"):
    print("  {:<22}: {}".format(
        name.lower() + " elements", len(unit_cell.element_sets[name])))
print("  boundary elements     : {}".format(
    len(unit_cell.element_sets["BOUNDARY_ELEMENTS"])))
print("  secondary elements    : {}".format(
    len(unit_cell.element_sets["SECONDARY_ELEMENTS"])))
print("  inner elements        : {}".format(
    len(unit_cell.element_sets["INNER_ELEMENTS"])))

mesh_file = os.path.join(here, "example_1.msh")
gmsh.write(mesh_file)
print("Periodic Gmsh mesh written to {}.".format(mesh_file))

# ---------------------------------------------------------------
# Transfer the reduced mixed-dimensional model to Abaqus. The explicit type
# map keeps the element choices visible in the example.
workdir = os.path.join(here, "abq")
if not os.path.isdir(workdir):
    os.makedirs(workdir)
os.chdir(workdir)

gmsh_unit_cell = unit_cell
element_types = {1: "B31", 3: "S4", 5: "C3D8"}
constituent_sets = ["CONTINUUM", "SHELL_X", "SHELL_Y", "SHELL_Z",
                    "BEAM_X", "BEAM_Y", "BEAM_Z"]

# Retain all periodic copies in a separate dummy model for visualization.
# This unreduced model is not used for homogenization.
periodic_representation = AbaqusRVE(
    gmsh_unit_cell, custom_types=element_types, remove_secondary=False)
assign_constituent_sections(periodic_representation, properties,
                            shell_thickness, beam_radius)
visualization_step = "VISUALIZATION"
periodic_representation.model.StaticStep(
    name=visualization_step, previous="Initial")


def zero_strain_history(time):
    return np.zeros((3, 3))


zero_history = LoadHistory(np.array([0.0, 1.0]), zero_strain_history)
periodic_representation.apply_periodic_BC(visualization_step, zero_history)
periodic_representation.run_analysis()
periodic_representation.export_to_vtk(
    os.path.join(here, "vtk_output", "periodic_representation"),
    field_outputs=["U"],
    element_sets=constituent_sets + ["SECONDARY_ELEMENTS"], clear=True)

# Rebuild the model with the ownership-consistent element reduction for the
# actual homogenization analyses.
unit_cell = AbaqusRVE(gmsh_unit_cell, custom_types=element_types)

# Check the ownership-consistent reduction before assigning sections.
elements_per_edge = nodes_per_edge - 1
expected_counts = {
    "CONTINUUM": elements_per_edge**3,
    "SHELL_X": elements_per_edge**2,
    "SHELL_Y": elements_per_edge**2,
    "SHELL_Z": elements_per_edge**2,
    "BEAM_X": elements_per_edge,
    "BEAM_Y": elements_per_edge,
    "BEAM_Z": elements_per_edge,
}

print("\nReduced Abaqus model:")
for name, expected in expected_counts.items():
    set_name = "ELSET_" + name
    actual = len(unit_cell.part.sets[set_name].elements)
    print("  {:<22}: {}".format(name.lower() + " elements", actual))
    if actual != expected:
        raise RuntimeError(
            "Element set {} contains {} elements; expected {}."
            .format(name, actual, expected))

expected_total = int(np.sum(list(expected_counts.values())))
actual_total = len(unit_cell.part.elements)
print("  total elements         : {}".format(actual_total))
if actual_total != expected_total:
    raise RuntimeError(
        "The reduced model contains {} elements; expected {}."
        .format(actual_total, expected_total))

assign_constituent_sections(unit_cell, properties, shell_thickness,
                            beam_radius)

step_name = "LOAD"
unit_cell.model.StaticStep(name=step_name, previous="Initial")

# ---------------------------------------------------------------
# Six independent strain cases determine the complete effective stiffness.
# The exact solution of the three normal-strain cases is affine, and their
# effective 3x3 normal-stiffness block is available in closed form. The shear
# cases additionally exercise the shell and beam rotational degrees of
# freedom in the Hill--Mandel check.
strain_magnitude = 0.01
volume = unit_cell.get_volume()
stiffness = np.zeros((6, 6))
hill_mandel_residuals = np.zeros(6)
energy_residuals = np.zeros(6)
power_tolerance = 1.0e-8
component_names = ("xx", "yy", "zz", "yz", "xz", "xy")

for component, component_name in enumerate(component_names):
    print("\n------------------------------")
    print("Run strain case {} ({}).".format(component + 1, component_name))

    def strain_history(time, index=component):
        strain = np.zeros(6)
        strain[index] = time*strain_magnitude
        return voigt_to_tensor(strain, "STRAIN")

    history = LoadHistory(np.array([0.0, 1.0]), strain_history)
    unit_cell.apply_periodic_BC(step_name, history)
    unit_cell.run_analysis(history_outputs=["ALLIE"])

    stress = unit_cell.average_stress(step_name)
    stress_voigt = tensor_to_voigt(stress, "STRESS")
    stiffness[:, component] = stress_voigt/strain_magnitude

    current_time = unit_cell.get_frame_time(step_name, -1)
    previous_time = unit_cell.get_frame_time(step_name, -2)
    strain_increment = (strain_history(current_time)
                        - strain_history(previous_time))
    macroscopic_power = float(np.sum(stress*strain_increment))
    microscopic_power = unit_cell.average_virtual_power(step_name)
    hill_mandel_residuals[component] = hill_mandel_residual(
        macroscopic_power, microscopic_power, power_tolerance)

    internal_energy_density = (_get_internal_energy(unit_cell, step_name)
                               / volume)
    macroscopic_energy_density = 0.5*float(np.sum(
        stress*strain_history(current_time)))
    energy_residuals[component] = _relative_residual(
        internal_energy_density, macroscopic_energy_density)

analytical_normal = analytical_normal_stiffness(
    edge_length, properties["CONTINUUM"], properties["SHELL"],
    properties["BEAM"], shell_thickness, beam_radius)
computed_normal = stiffness[:3, :3]
relative_error = (np.linalg.norm(computed_normal - analytical_normal)
                  / np.linalg.norm(analytical_normal))
major_symmetry_residual = (np.max(np.abs(stiffness - stiffness.T))
                           / np.max(np.abs(stiffness)))

print("\n==============================")
print("Computed effective stiffness:")
print(np.array2string(stiffness, precision=8, suppress_small=True))
print("\nAnalytical normal-stiffness block:")
print(np.array2string(analytical_normal, precision=8, suppress_small=True))
print("\nRelative normal-block error : {:.3e}".format(relative_error))
print("Major-symmetry residual     : {:.3e}".format(
    major_symmetry_residual))
print("Cubic-symmetry residual     : {:.3e}".format(
    cubic_symmetry_residual(stiffness)))
print("Hill--Mandel residuals      : {}".format(
    np.array2string(hill_mandel_residuals, precision=3)))
print("Energy residuals            : {}".format(
    np.array2string(energy_residuals, precision=3)))

result_dir = os.path.join(here, "results")
if not os.path.isdir(result_dir):
    os.makedirs(result_dir)
export_to_csv_file(stiffness,
                   os.path.join(result_dir, "effective_stiffness_computed"))
export_to_csv_file(computed_normal,
                   os.path.join(result_dir, "normal_stiffness_computed"))
export_to_csv_file(analytical_normal,
                   os.path.join(result_dir, "normal_stiffness_analytical"))
diagnostics = np.column_stack((hill_mandel_residuals, energy_residuals))
export_to_csv_file(diagnostics,
                   os.path.join(result_dir, "consistency_residuals"))

# The final ODB contains the xy-shear case and provides the deformed reduced
# representation for visualization alongside the complete periodic model.
unit_cell.export_to_vtk(
    os.path.join(here, "vtk_output", "reduced_xy_shear"),
    field_outputs=["U", "UR"], element_sets=constituent_sets, clear=True)

# ---------------------------------------------------------------
# Independently verify the reduction against the complete periodic cell. For
# the present linear-elastic materials, scaling the shell and beam moduli by
# their ownership factors is equivalent to weighting their complete element
# virtual work by one half and one quarter, respectively.
weighted_workdir = os.path.join(here, "abq_weighted")
if not os.path.isdir(weighted_workdir):
    os.makedirs(weighted_workdir)
os.chdir(weighted_workdir)

weighted_cell = AbaqusRVE(
    gmsh_unit_cell, custom_types=element_types, remove_secondary=False)
gmsh.finalize()

expected_full_counts = {
    "CONTINUUM": elements_per_edge**3,
    "SHELL_X": 2*elements_per_edge**2,
    "SHELL_Y": 2*elements_per_edge**2,
    "SHELL_Z": 2*elements_per_edge**2,
    "BEAM_X": 4*elements_per_edge,
    "BEAM_Y": 4*elements_per_edge,
    "BEAM_Z": 4*elements_per_edge,
}

print("\nComplete ownership-weighted Abaqus model:")
for name, expected in expected_full_counts.items():
    set_name = "ELSET_" + name
    actual = len(weighted_cell.part.sets[set_name].elements)
    print("  {:<22}: {}".format(name.lower() + " elements", actual))
    if actual != expected:
        raise RuntimeError(
            "Element set {} contains {} elements; expected {}."
            .format(name, actual, expected))

expected_full_total = int(np.sum(list(expected_full_counts.values())))
actual_full_total = len(weighted_cell.part.elements)
print("  total elements         : {}".format(actual_full_total))
if actual_full_total != expected_full_total:
    raise RuntimeError(
        "The complete model contains {} elements; expected {}."
        .format(actual_full_total, expected_full_total))

weighted_properties = {
    "CONTINUUM": properties["CONTINUUM"],
    "SHELL": (0.5*properties["SHELL"][0], properties["SHELL"][1]),
    "BEAM": (0.25*properties["BEAM"][0], properties["BEAM"][1]),
}
assign_constituent_sections(weighted_cell, weighted_properties,
                            shell_thickness, beam_radius)

weighted_cell.model.StaticStep(name=step_name, previous="Initial")
weighted_stiffness = np.zeros((6, 6))

for component, component_name in enumerate(component_names):
    print("\n------------------------------")
    print("Run weighted strain case {} ({}).".format(
        component + 1, component_name))

    def weighted_strain_history(time, index=component):
        strain = np.zeros(6)
        strain[index] = time*strain_magnitude
        return voigt_to_tensor(strain, "STRAIN")

    weighted_history = LoadHistory(
        np.array([0.0, 1.0]), weighted_strain_history)
    weighted_cell.apply_periodic_BC(step_name, weighted_history)
    weighted_cell.run_analysis()

    weighted_stress = weighted_cell.average_stress(step_name)
    weighted_stiffness[:, component] = (
        tensor_to_voigt(weighted_stress, "STRESS")/strain_magnitude)

reduction_residual = (
    np.linalg.norm(stiffness - weighted_stiffness)
    / np.linalg.norm(stiffness))

print("\n==============================")
print("Complete ownership-weighted effective stiffness:")
print(np.array2string(weighted_stiffness, precision=8,
                      suppress_small=True))
print("\nReduction residual         : {:.3e}".format(
    reduction_residual))

export_to_csv_file(
    weighted_stiffness,
    os.path.join(result_dir, "effective_stiffness_weighted"))
export_to_csv_file(
    np.array([[reduction_residual]]),
    os.path.join(result_dir, "reduction_residual"))
