"""
Generate a periodic sheet-gyroid surface mesh with microgen and export it
to 'gyroid.msh' for the Abaqus example (main.py).

Run this once, separately, in a Python environment that has microgen
installed (Abaqus itself does not need microgen -- main.py only reads the
resulting mesh file):

    python generate_mesh.py

microgen builds the gyroid mid-surface (phi = 0) as a triangulated pyvista
mesh. The cell is phase-shifted slightly so its corners avoid the surface,
which makes the mesh exactly periodic across opposite faces (checked with
microgen.is_periodic). The surface is written as a Gmsh discrete surface
with the physical group 'GYROID', which main.py imports via
GmshRVE.read_from_mesh.
"""

import os
import inspect

import numpy as np

# microgen 1.2.1 still uses a few numpy aliases removed in numpy 2.0
# (e.g. np.float_). Shim them so it imports under a numpy>=2 environment.
for _old, _new in [("float_", "float64"), ("complex_", "complex128"),
                   ("unicode_", "str_"), ("bool8", "bool_")]:
    if not hasattr(np, _old):
        try:
            setattr(np, _old, getattr(np, _new))
        except AttributeError:
            pass

import gmsh
import microgen
from microgen.shape.surface_functions import gyroid

here = os.path.dirname(os.path.abspath(
    inspect.getframeinfo(inspect.currentframe()).filename))

# Geometry parameters
cell_size = 1.0      # edge length of the cubic unit cell
resolution = 30      # grid resolution (finer meshes need a full Abaqus
                     # license; reduce it for the node-limited Learning Edition)
phase_shift = 0.3    # shift so the cube corners avoid the surface

# Build the gyroid mid-surface (phi = 0) with microgen
tpms = microgen.Tpms(surface_function=gyroid, offset=0.0, cell_size=cell_size,
                     resolution=resolution,
                     phase_shift=(phase_shift, phase_shift, phase_shift))
surface = tpms.surface
points = np.asarray(surface.points)
if not microgen.is_periodic(points):
    raise RuntimeError(
        "The generated microgen surface is not periodic; adjust phase_shift.")

# pyvista stores an all-triangle surface as [3, i, j, k, 3, i, j, k, ...]
triangles = surface.faces.reshape(-1, 4)[:, 1:]

# Write the surface as a Gmsh discrete surface with a physical group
gmsh.initialize()
gmsh.model.add("GYROID")
entity = gmsh.model.add_discrete_entity(2)
gmsh.model.mesh.add_nodes(2, entity, list(range(1, len(points) + 1)),
                          points.flatten().tolist())
gmsh.model.mesh.add_elements_by_type(entity, 2,
                                     list(range(1, len(triangles) + 1)),
                                     (triangles + 1).flatten().tolist())
group = gmsh.model.add_physical_group(2, [entity])
gmsh.model.set_physical_name(2, group, "GYROID")

mesh_file = os.path.join(here, "gyroid.msh")
gmsh.option.set_number("Mesh.MshFileVersion", 4.1)
gmsh.write(mesh_file)
gmsh.finalize()

print("Wrote {} ({} nodes, {} triangles).".format(
    mesh_file, len(points), len(triangles)))
