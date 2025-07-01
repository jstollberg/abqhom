import os
import paraview.simple as simple


def find_all_pvd_files(directory):
    pvd_files = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".pvd"):
                pvd_files.append(os.path.join(root, file))

    return pvd_files


def visualize_deformed_RVE(pvd_file):
    # Load PVD file
    pvd = simple.OpenDataFile(pvd_file)
    pvd.UpdatePipeline()
    info = pvd.GetDataInformation()

    # Check if the model contains a deformation field
    deformation = False
    point_data = info.GetPointDataInformation()
    for i in range(point_data.GetNumberOfArrays()):
        field_name = point_data.GetArrayInformation(i).GetName()
        if field_name == "U":
            deformation = True

    # Get cell field names and check for different physical groups
    radius = False
    group_names = []
    cell_data = info.GetCellDataInformation()
    for i in range(cell_data.GetNumberOfArrays()):
        field_name = cell_data.GetArrayInformation(i).GetName()

        # Check if the model contains radius information
        if field_name == "RADIUS":
            radius = True

        # Collect physical groups
        if "_group" in field_name:
            group_names.append(field_name)

    # Loop over all physical groups
    displays = []
    groups = []
    for threshold_field in group_names:

        # Apply threshold for the current physical group
        threshold = simple.Threshold(Input=pvd)
        threshold.UpdatePipeline()
        threshold.Scalars = ["CELLS", threshold_field]
        threshold.LowerThreshold = 1.0
        threshold.UpperThreshold = 1.0

        groups.append(threshold)

    # Connect the physical groups
    group_data_sets = simple.GroupDatasets(Input=groups)
    group_data_sets.UpdatePipeline()

    # Visualize beams or trusses as tubes
    tube = group_data_sets
    if radius:
        cell_to_point = simple.CellDatatoPointData(Input=group_data_sets)
        cell_to_point.UpdatePipeline()

        surface = simple.ExtractSurface(Input=cell_to_point)
        surface.UpdatePipeline()
        tube = simple.Tube(Input=surface)
        tube.UpdatePipeline()
        tube.Scalars = ["POINTS", "RADIUS"]
        tube.VaryRadius = "By Absolute Scalar"
        tube.NumberofSides = 50

    # Show the deformed structure if there is a displacement field
    warp = tube
    if deformation:
        warp = simple.WarpByVector(Input=tube)
        warp.UpdatePipeline()
        warp.Vectors = ["POINTS", "U"]

    # Show the last time step (usually the deformed state)
    animationScene = simple.GetAnimationScene()
    animationScene.UpdateAnimationUsingDataTimeSteps()
    timesteps = animationScene.TimeKeeper.TimestepValues
    animationScene.AnimationTime = timesteps[-1]
    warp.UpdatePipeline(time=timesteps[-1])

    # Show results
    display = simple.Show(warp)
    if deformation:
        simple.ColorBy(display, ["POINTS", "U"])
        display.RescaleTransferFunctionToDataRange(True)
        display.SetScalarBarVisibility(simple.GetActiveView(), True)

    displays.append(display)

    simple.Render()
    return displays


if __name__ == "__main__":
    # ---------------------------------------
    # Define the folder containing the files
    # ---------------------------------------
    vtk_folder = None
    # ---------------------------------------

    # Find all pvd files
    if vtk_folder is None:
        vtk_folder = os.getcwd()
    pvd_files = find_all_pvd_files(vtk_folder)

    # Visualize all pvd files
    all_displays = []
    for file in pvd_files:
        all_displays.append(visualize_deformed_RVE(file))

    # Hide all displays except the last one
    for display_group in all_displays[0:-1]:
        for display in display_group:
            display.Visibility = 0
