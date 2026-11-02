# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-or-later

"""Generate the periodic sheet-gyroid mesh used by ``main.py``.

Microgen constructs the gyroid midsurface as a triangulated PyVista mesh. A
phase shift keeps the surface away from the cell corners, and periodicity is
verified before the mesh is written as a Gmsh discrete surface in
``gyroid.msh``.
"""

import os

import numpy as np

import gmsh
import microgen
from microgen.shape.surface_functions import gyroid

def main():
    """Generate and write the periodic sheet-gyroid mesh."""
    try:
        here = os.getcwd()

        # Geometry parameters
        cell_size = 1.0      # edge length of the cubic unit cell
        resolution = 30      # grid resolution
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

        print("Wrote {} ({} nodes, {} triangles).".format(
            mesh_file, len(points), len(triangles)))
    finally:
        if gmsh.is_initialized():
            gmsh.finalize()


if __name__ == "__main__":
    main()
