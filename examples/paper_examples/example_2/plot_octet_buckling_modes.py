"""Plot the first two periodic buckling modes of the Octet cell."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import vtk
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from vtk.util.numpy_support import vtk_to_numpy


ROOT = Path(__file__).resolve().parent
VTU_DIRECTORY = (
    ROOT
    / "final_vtk"
    / "buckling"
    / "EXAMPLE_2_OCTET_BUCKLING"
    / "vtu"
)
RESULT_DIRECTORY = ROOT / "final_results" / "buckling"
DISPLAY_SCALE = 0.12


def read_mode(mode):
    path = VTU_DIRECTORY / f"EXAMPLE_2_OCTET_BUCKLING_0_{mode}.vtu"
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = reader.GetOutput()

    points = vtk_to_numpy(grid.GetPoints().GetData())
    displacement = vtk_to_numpy(grid.GetPointData().GetArray("U"))
    connectivity = []
    for cell_index in range(grid.GetNumberOfCells()):
        cell = grid.GetCell(cell_index)
        connectivity.append(
            [cell.GetPointId(local_index) for local_index in range(cell.GetNumberOfPoints())]
        )
    return points, displacement, connectivity


def segments(points, connectivity):
    return np.asarray(
        [
            [points[cell[local_index]], points[cell[local_index + 1]]]
            for cell in connectivity
            for local_index in range(len(cell) - 1)
        ]
    )


def set_equal_axes(axis):
    axis.set_xlim(-0.65, 0.65)
    axis.set_ylim(-0.65, 0.65)
    axis.set_zlim(-0.65, 0.65)
    axis.set_box_aspect((1.0, 1.0, 1.0))
    axis.set_xlabel(r"$X_1$")
    axis.set_ylabel(r"$X_2$")
    axis.set_zlabel(r"$X_3$")
    axis.view_init(elev=22.0, azim=-55.0)


def main():
    eigenvalues = np.loadtxt(RESULT_DIRECTORY / "eigenvalues.csv", delimiter=",")
    figure = plt.figure(figsize=(8.0, 3.7), constrained_layout=True)

    color_map = plt.get_cmap("viridis")
    for panel, mode in enumerate((1, 2), start=1):
        points, displacement, connectivity = read_mode(mode)
        deformed = points + DISPLAY_SCALE * displacement
        undeformed_segments = segments(points, connectivity)
        deformed_segments = segments(deformed, connectivity)

        axis = figure.add_subplot(1, 2, panel, projection="3d")
        axis.add_collection3d(
            Line3DCollection(
                undeformed_segments,
                colors="0.78",
                linewidths=0.55,
                alpha=0.75,
            )
        )

        magnitude = np.linalg.norm(displacement, axis=1)
        segment_values = np.asarray(
            [
                0.5 * (magnitude[cell[i]] + magnitude[cell[i + 1]])
                for cell in connectivity
                for i in range(len(cell) - 1)
            ]
        )
        collection = Line3DCollection(
            deformed_segments,
            cmap=color_map,
            linewidths=1.8,
        )
        collection.set_array(segment_values)
        collection.set_clim(0.0, magnitude.max())
        axis.add_collection3d(collection)
        set_equal_axes(axis)
        eigenvalue = float(eigenvalues[mode - 1, 1])
        axis.set_title(
            rf"Mode {mode}: $\lambda={eigenvalue:.5f}$",
            fontsize=10,
        )

    output = RESULT_DIRECTORY / "first_buckling_modes.png"
    figure.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(figure)
    print(output)


if __name__ == "__main__":
    main()
