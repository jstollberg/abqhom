import os
import sys
import numpy as np

base_path = "C:/Users/jonat/Documents/"
gmsh_path = os.path.join(base_path, "gmsh", "lib")
abqhom_path = os.path.join(base_path, "abqhom", "src")
paraqus_path = os.path.join(base_path, "Paraqus", "src")
for i in (gmsh_path, abqhom_path, paraqus_path):
    sys.path.append(i)

# Import Gmsh to create the geometry
import gmsh

# Import homogenization tools
from abqhom.utils import voigt_to_tensor, tensor_to_voigt
from abqhom.gmshmodel import GmshRVE
from abqhom.abqmodel import AbaqusRVE, LoadHistory
from abqhom.utils import reuss, voigt, export_to_csv_file


def create_geometry(model_name, edge_length=1.0, nodes_matrix=10,
                    nodes_fibers=3):

    # Initialize Gmsh and add a new model
    if not gmsh.is_initialized():
        gmsh.initialize()
    gmsh.model.add(model_name)

    gmsh.option.set_number("Geometry.Tolerance", 1e-5)

    x = y = z = 0.0  # Center point
    dx = dy = dz = edge_length  # Edge length

    # Add a cube to the model representing the matrix
    p1 = gmsh.model.geo.add_point(x - 0.5*dx, y - 0.5*dy, z - 0.5*dz)
    p2 = gmsh.model.geo.add_point(x - 0.5*dx, y - 0.5*dy, z + 0.5*dz)
    p3 = gmsh.model.geo.add_point(x + 0.5*dx, y - 0.5*dy, z + 0.5*dz)
    p4 = gmsh.model.geo.add_point(x + 0.5*dx, y - 0.5*dy, z - 0.5*dz)
    p5 = gmsh.model.geo.add_point(x - 0.5*dx, y + 0.5*dy, z - 0.5*dz)
    p6 = gmsh.model.geo.add_point(x - 0.5*dx, y + 0.5*dy, z + 0.5*dz)
    p7 = gmsh.model.geo.add_point(x + 0.5*dx, y + 0.5*dy, z + 0.5*dz)
    p8 = gmsh.model.geo.add_point(x + 0.5*dx, y + 0.5*dy, z - 0.5*dz)

    line_connectivity = [[p1, p2], [p2, p3], [p3, p4], [p4, p1],
                         [p5, p6], [p6, p7], [p7, p8], [p8, p5],
                         [p1, p5], [p2, p6], [p3, p7], [p4, p8]]
    lines = [gmsh.model.geo.add_line(*i) for i in line_connectivity]

    surfaces = []
    surface_connectivity = [[lines[0], lines[1], lines[2], lines[3]],
                            [lines[4], lines[5], lines[6], lines[7]],
                            [lines[0], lines[9], lines[4], lines[8]],
                            [lines[1], lines[10], lines[5], lines[9]],
                            [lines[2], lines[11], lines[6], lines[10]],
                            [lines[3], lines[8], lines[7], lines[11]]]
    for i in surface_connectivity:
        loop = gmsh.model.geo.add_curve_loop(i, reorient=True)
        surfaces.append(gmsh.model.geo.add_plane_surface([loop]))

    loop = gmsh.model.geo.add_surface_loop(surfaces)
    volume = gmsh.model.geo.add_volume([loop])

    # Add transfinite hex mesh constraints on the cube
    for i in lines:
        gmsh.model.geo.mesh.set_transfinite_curve(i, nodes_matrix)

    for i in surfaces:
        gmsh.model.geo.mesh.set_transfinite_surface(i)
        gmsh.model.geo.mesh.set_recombine(2, i)

    gmsh.model.geo.mesh.set_transfinite_volume(volume)

    # Add one fiber
    # TODO: Here you need to generate your fibers
    fiber_p1 = gmsh.model.geo.add_point(-0.25*dx, 0.25*dy, -0.5*dz + 0.001*dz)
    fiber_p2 = gmsh.model.geo.add_point(-0.25*dx, 0.25*dy, 0.5*dz - 0.001*dz)
    fiber_p3 = gmsh.model.geo.add_point(0.25*dx, -0.25*dy, -0.5*dz + 0.001*dz)
    fiber_p4 = gmsh.model.geo.add_point(0.25*dx, -0.25*dy, 0.5*dz - 0.001*dz)
    fiber_1 = gmsh.model.geo.add_line(fiber_p1, fiber_p2)
    fiber_2 = gmsh.model.geo.add_line(fiber_p3, fiber_p4)
    gmsh.model.geo.mesh.set_transfinite_curve(fiber_1, nodes_fibers)
    gmsh.model.geo.mesh.set_transfinite_curve(fiber_2, nodes_fibers)
    fibers = [fiber_1, fiber_2]

    # Synchronize the geo kernel to save the geometry
    gmsh.model.geo.synchronize()

    # Add physical groups
    matrix = gmsh.model.add_physical_group(3, [volume])
    gmsh.model.set_physical_name(3, matrix, "MATRIX")
    fibers = gmsh.model.add_physical_group(1, fibers)
    gmsh.model.set_physical_name(1, fibers, "FIBER")


# ---------------------------------------------------------------
# Set working directory
workdir = os.path.join(base_path,
                       "abqhom/examples/3d_1d_coupled/abq")
if not os.path.isdir(workdir):
    os.mkdir(workdir)
os.chdir(workdir)

# Name folder to store VTK files
vtk_path = os.path.join(workdir, "..", "vtk_output")

# Material parameters: Young's modulus - Poisson's ratio
parameters = {"MATRIX": [3.5, 0.35],
              "FIBER":  [230.0, 0.2]}

# ---------------------------------------------------------------
# Run homogenization routine
C = np.empty((6, 6))
for i in range(6):
    print("\n------------------------------")
    print("Run case {}.".format(i + 1))

    # Define load history
    strain_magnitude = 0.1 if i < 3 else 0.2

    def strain_history(time):
        strain = np.zeros(6)
        strain[i] = strain_magnitude
        strain = voigt_to_tensor(strain, "STRAIN")
        return time*strain

    history = LoadHistory(np.array([0.0, 1.0]), strain_history)
    print("Strain = " + str(tensor_to_voigt(history(1.0), "STRAIN")
                            .round(3).tolist()))

    # Create geometry and periodic mesh in Gmsh
    model_name = "RVE_{}".format(i + 1)
    create_geometry(model_name, nodes_matrix=9, nodes_fibers=20)
    unit_cell = GmshRVE(model_name)
    unit_cell.add_periodic_mesh_constraints()
    unit_cell.generate_mesh()

    # Load the model into Abaqus
    unit_cell = AbaqusRVE(unit_cell)

    # Assign material
    for name in parameters:
        E, nu = parameters[name]
        material = unit_cell.model.Material(name=name)
        material.Elastic(table=[[E, nu]])
    unit_cell.assign_solid_section("MATRIX", "MATRIX")
    unit_cell.assign_beam_section("FIBER", "FIBER", 0.02)

    # Activate embedded region constraint for fibers
    matrix_elements = unit_cell.instance.sets["ELSET_MATRIX"]
    fiber_elements = unit_cell.instance.sets["ELSET_FIBER"]
    unit_cell.model.EmbeddedRegion("COUPLING", fiber_elements, matrix_elements)

    # Add load step and boundary conditions
    unit_cell.model.StaticStep(name="LOAD", previous="Initial")
    unit_cell.apply_periodic_BC("LOAD", history)

    # Run FE analysis
    unit_cell.run_analysis(field_outputs=["NFORC", "U", "S"])

    # Export to VTK
    unit_cell.export_to_vtk(vtk_path, field_outputs=["U", ("S", "MISES")],
                            element_sets=["MATRIX", "FIBER"])

    gmsh.finalize()
