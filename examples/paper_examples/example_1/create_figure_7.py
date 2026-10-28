"""Create Figure 7: complete and reduced mixed-dimensional unit cells.

VTK renders the finite element models, while Matplotlib supplies the panel
labels and categorical legend so that all typography matches the manuscript.
Run the script with the Python interpreter of the ``chatgpt`` environment.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.transforms import Bbox
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk


HERE = Path(__file__).resolve().parent
REPOSITORY = HERE.parents[3]
VTK_ROOT = HERE / "vtk_output"
OUTPUT_PDF = REPOSITORY / "tex_source_elsevier" / "figure_7.pdf"
OUTPUT_PNG = HERE / "figure_7.png"

COMPLETE_FILE = (
    VTK_ROOT
    / "periodic_representation"
    / "EXAMPLE_1"
    / "vtu"
    / "EXAMPLE_1_0_1.vtu"
)
REDUCED_FILE = (
    VTK_ROOT
    / "reduced_xy_shear"
    / "EXAMPLE_1"
    / "vtu"
    / "EXAMPLE_1_0_0.vtu"
)

# Constituent colors matching Figure 3 of the manuscript (continuum = silver
# matrix of Fig. 3c, shell = teal, beam = slate).
COLORS = {
    "continuum": "#c7ccd4",
    "shell": "#4f9a8e",
    "beam": "#627489",
}

# S4 shell thickness (Table 1); the flat QUAD cells are extruded to a slab of
# this thickness so the shells render with their physical thickness.
SHELL_THICKNESS = 0.02

# elsarticle 5p/twocolumn metrics. matplotlib pt are 1/72 in, LaTeX pt are
# 1/72.27 in; authoring the figure at its on-page display width makes a label
# set at FS10 render as a true 10 pt once included at 0.95\columnwidth.
COLUMN_WIDTH_IN = 252.0 / 72.27
FS8 = 8.0 * 72.0 / 72.27
FS10 = 10.0 * 72.0 / 72.27


def hex_to_rgb(color):
    """Return an RGB tuple on the interval [0, 1]."""
    color = color.lstrip("#")
    return tuple(int(color[i:i + 2], 16) / 255.0 for i in (0, 2, 4))


def read_grid(path):
    """Read an unstructured VTK grid."""
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = vtk.vtkUnstructuredGrid()
    grid.ShallowCopy(reader.GetOutput())
    return grid


def extract_cell_type(grid, cell_type):
    """Extract all cells of one VTK cell type."""
    ids = vtk.vtkIdList()
    for cell_id in range(grid.GetNumberOfCells()):
        if grid.GetCellType(cell_id) == cell_type:
            ids.InsertNextId(cell_id)

    extractor = vtk.vtkExtractCells()
    extractor.SetInputData(grid)
    extractor.SetCellList(ids)
    extractor.Update()

    extracted = vtk.vtkUnstructuredGrid()
    extracted.ShallowCopy(extractor.GetOutput())
    return extracted


def radially_scale(grid, factor):
    """Apply a minute radial offset to avoid coincident-surface artifacts."""
    scaled = vtk.vtkUnstructuredGrid()
    scaled.DeepCopy(grid)
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(scaled.GetNumberOfPoints())
    for point_id in range(scaled.GetNumberOfPoints()):
        point = np.asarray(scaled.GetPoint(point_id)) * factor
        points.SetPoint(point_id, point)
    scaled.SetPoints(points)
    return scaled


def surface_actor(grid, face_color, opacity, specular, ambient, diffuse):
    """Smooth-shaded opaque surface actor (Figure 3 style, no element edges)."""
    surface = vtk.vtkDataSetSurfaceFilter()
    surface.SetInputData(grid)

    normals = vtk.vtkPolyDataNormals()
    normals.SetInputConnection(surface.GetOutputPort())
    normals.SetFeatureAngle(60.0)
    normals.SplittingOff()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(normals.GetOutputPort())
    mapper.ScalarVisibilityOff()

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    prop = actor.GetProperty()
    prop.SetColor(hex_to_rgb(face_color))
    prop.SetOpacity(opacity)
    prop.SetInterpolationToPhong()
    prop.SetAmbient(ambient)
    prop.SetDiffuse(diffuse)
    prop.SetSpecular(specular)
    prop.SetSpecularPower(20.0)
    prop.EdgeVisibilityOff()
    return actor


def extrude_shell_surface(grid, thickness, center):
    """Extrude the flat shell quads outward into a slab of the given thickness.

    Point normals are split at the cube edges and oriented away from the cell
    center, so each face extrudes along its own outward normal.
    """
    surface = vtk.vtkDataSetSurfaceFilter()
    surface.SetInputData(grid)

    normals = vtk.vtkPolyDataNormals()
    normals.SetInputConnection(surface.GetOutputPort())
    normals.SplittingOn()
    normals.SetFeatureAngle(30.0)
    normals.ConsistencyOn()
    normals.ComputePointNormalsOn()
    normals.ComputeCellNormalsOff()
    normals.Update()

    poly = vtk.vtkPolyData()
    poly.DeepCopy(normals.GetOutput())
    points = vtk_to_numpy(poly.GetPoints().GetData())
    point_normals = vtk_to_numpy(poly.GetPointData().GetNormals()).copy()
    outward = points - np.asarray(center, dtype=float)
    flip = np.einsum("ij,ij->i", point_normals, outward) < 0.0
    point_normals[flip] = -point_normals[flip]
    poly.GetPointData().SetNormals(
        numpy_to_vtk(np.ascontiguousarray(point_normals), deep=True))

    extrude = vtk.vtkLinearExtrusionFilter()
    extrude.SetInputData(poly)
    extrude.SetExtrusionTypeToNormalExtrusion()
    extrude.SetScaleFactor(thickness)
    extrude.CappingOn()

    out_normals = vtk.vtkPolyDataNormals()
    out_normals.SetInputConnection(extrude.GetOutputPort())
    out_normals.SetFeatureAngle(45.0)
    out_normals.SplittingOn()
    out_normals.ConsistencyOn()
    out_normals.Update()
    return out_normals.GetOutput()


def shell_slab_actor(grid, face_color, thickness, center,
                     specular, ambient, diffuse):
    """Smooth-shaded opaque shell slab of the physical shell thickness."""
    slab = extrude_shell_surface(grid, thickness, center)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(slab)
    mapper.ScalarVisibilityOff()

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    prop = actor.GetProperty()
    prop.SetColor(hex_to_rgb(face_color))
    prop.SetInterpolationToPhong()
    prop.SetAmbient(ambient)
    prop.SetDiffuse(diffuse)
    prop.SetSpecular(specular)
    prop.SetSpecularPower(20.0)
    prop.EdgeVisibilityOff()
    return actor


def beam_actor(grid):
    """Render beam centerlines as tubes using the actual section radius."""
    radius = 0.04
    radius_array = grid.GetCellData().GetArray("RADIUS")
    if radius_array is not None and radius_array.GetNumberOfTuples() > 0:
        radius = float(radius_array.GetRange()[1])

    geometry = vtk.vtkGeometryFilter()
    geometry.SetInputData(grid)

    tube = vtk.vtkTubeFilter()
    tube.SetInputConnection(geometry.GetOutputPort())
    tube.SetRadius(radius)
    tube.SetNumberOfSides(24)
    tube.CappingOn()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(tube.GetOutputPort())
    mapper.ScalarVisibilityOff()

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(hex_to_rgb(COLORS["beam"]))
    actor.GetProperty().SetInterpolationToPhong()
    actor.GetProperty().SetAmbient(0.30)
    actor.GetProperty().SetDiffuse(0.75)
    actor.GetProperty().SetSpecular(0.40)
    actor.GetProperty().SetSpecularPower(20.0)
    return actor


def render_cell(path, image_size=(1500, 1250)):
    """Render one cell from a common isometric camera."""
    grid = read_grid(path)
    continuum = extract_cell_type(grid, vtk.VTK_HEXAHEDRON)
    shells = radially_scale(extract_cell_type(grid, vtk.VTK_QUAD), 1.004)
    beams = radially_scale(extract_cell_type(grid, vtk.VTK_LINE), 1.010)

    bounds = continuum.GetBounds()
    center = (0.5 * (bounds[0] + bounds[1]),
              0.5 * (bounds[2] + bounds[3]),
              0.5 * (bounds[4] + bounds[5]))

    renderer = vtk.vtkRenderer()
    renderer.SetBackground(1.0, 1.0, 1.0)
    renderer.SetBackgroundAlpha(0.0)

    # Three-light setup for the smooth shaded look of Figure 3.
    light_kit = vtk.vtkLightKit()
    light_kit.AddLightsToRenderer(renderer)

    renderer.AddActor(surface_actor(
        continuum, COLORS["continuum"],
        opacity=1.0, specular=0.15, ambient=0.35, diffuse=0.75))
    renderer.AddActor(shell_slab_actor(
        shells, COLORS["shell"], SHELL_THICKNESS, center,
        specular=0.30, ambient=0.30, diffuse=0.80))
    renderer.AddActor(beam_actor(beams))

    camera = vtk.vtkCamera()
    camera.SetPosition(-2.35, -2.55, 2.00)
    camera.SetFocalPoint(0.0, 0.0, 0.0)
    camera.SetViewUp(0.0, 0.0, 1.0)
    camera.ParallelProjectionOn()
    camera.SetParallelScale(0.88)
    renderer.SetActiveCamera(camera)

    window = vtk.vtkRenderWindow()
    window.SetOffScreenRendering(True)
    window.SetAlphaBitPlanes(True)
    window.SetMultiSamples(0)
    window.SetSize(*image_size)
    window.AddRenderer(renderer)
    window.Render()

    capture = vtk.vtkWindowToImageFilter()
    capture.SetInput(window)
    capture.SetInputBufferTypeToRGBA()
    capture.ReadFrontBufferOff()
    capture.Update()

    image = capture.GetOutput()
    width, height, _ = image.GetDimensions()
    rgba = vtk_to_numpy(image.GetPointData().GetScalars())
    rgba = rgba.reshape(height, width, 4)
    rgba = np.flipud(rgba)
    window.Finalize()
    return rgba


def crop_to_content(rgba, pad=12, thr=8):
    """Trim transparent margins so the panels sit close together."""
    alpha = rgba[:, :, 3]
    rows = np.where(alpha.max(axis=1) > thr)[0]
    cols = np.where(alpha.max(axis=0) > thr)[0]
    if rows.size == 0 or cols.size == 0:
        return rgba
    r0 = max(rows[0] - pad, 0)
    r1 = min(rows[-1] + 1 + pad, rgba.shape[0])
    c0 = max(cols[0] - pad, 0)
    c1 = min(cols[-1] + 1 + pad, rgba.shape[1])
    return rgba[r0:r1, c0:c1]


def configure_matplotlib():
    """Use manuscript-compatible serif typography."""
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "Nimbus Roman",
                       "DejaVu Serif"],
        "font.size": FS10,
        "mathtext.fontset": "stix",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def create_figure():
    """Assemble and export the single-column two-panel publication figure.

    Authored at 0.95\\columnwidth with a fixed figure width (no tight bbox), so
    the embedded PDF keeps that exact natural width and all text renders at a
    true 10 pt when included at ``width=0.95\\columnwidth``.
    """
    configure_matplotlib()
    complete = crop_to_content(render_cell(COMPLETE_FILE))
    reduced = crop_to_content(render_cell(REDUCED_FILE))

    width_in = 0.95 * COLUMN_WIDTH_IN
    side = 0.010 * width_in
    center_gap = 0.040 * width_in
    w_cell = (width_in - 2.0 * side - center_gap) / 2.0

    def aspect(image):  # width / height in pixels
        return image.shape[1] / image.shape[0]

    images = (complete, reduced)
    h_cell = max(w_cell / aspect(image) for image in images)

    top_pad = 0.02
    legend_h = 0.24
    height_in = top_pad + h_cell + legend_h

    figure = plt.figure(figsize=(width_in, height_in))
    for i, (image, panel) in enumerate(zip(images, ("(a)", "(b)"))):
        h = w_cell / aspect(image)
        cx = side + i * (w_cell + center_gap) + w_cell / 2.0
        axis = figure.add_axes([
            (cx - w_cell / 2.0) / width_in,
            (legend_h + (h_cell - h)) / height_in,
            w_cell / width_in,
            h / height_in,
        ])
        axis.imshow(image)
        axis.set_axis_off()
        axis.text(0.01, 0.99, panel, transform=axis.transAxes,
                  ha="left", va="top", fontsize=FS10)

    handles = [
        Patch(facecolor=COLORS["continuum"], edgecolor="none", label="Continuum"),
        Patch(facecolor=COLORS["shell"], edgecolor="none", label="Shell"),
        Patch(facecolor=COLORS["beam"], edgecolor="none", label="Beam"),
    ]
    figure.legend(
        handles=handles,
        loc="lower center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 0.005),
        handlelength=1.2,
        handletextpad=0.5,
        columnspacing=1.6,
        borderpad=0.0,
        fontsize=FS10,
    )

    # Crop only the whitespace below the legend, keeping the full authored
    # width so the embedded PDF stays 0.95\columnwidth and the text true 8 pt.
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    legend = figure.legends[0]
    leg_bbox = legend.get_window_extent(renderer).transformed(
        figure.dpi_scale_trans.inverted())
    bottom = max(leg_bbox.y0 - top_pad, 0.0)
    crop = Bbox([[0.0, bottom], [width_in, height_in]])

    figure.savefig(OUTPUT_PDF, dpi=600, transparent=True, bbox_inches=crop)
    figure.savefig(OUTPUT_PNG, dpi=600, bbox_inches=crop)
    plt.close(figure)

    print("Wrote {}".format(OUTPUT_PDF))
    print("Wrote {}".format(OUTPUT_PNG))


if __name__ == "__main__":
    create_figure()
