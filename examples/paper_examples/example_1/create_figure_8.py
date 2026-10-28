"""Create Figure 8: reduced mixed-dimensional cell under the 12-shear state.

VTK renders the warped reduced cell colored by the magnitude of the nodal
rotation vector (Abaqus UR). Matplotlib supplies a horizontal colorbar below
the cell so that all typography matches the manuscript. The VTK lookup table
is sampled from the same Matplotlib colormap and normalization as the
colorbar, so the two agree exactly.

Run the script with the Python interpreter of the ``chatgpt`` environment.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from matplotlib.ticker import FormatStrFormatter
from matplotlib.transforms import Bbox
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk


HERE = Path(__file__).resolve().parent
REPOSITORY = HERE.parents[3]
VTK_ROOT = HERE / "vtk_output"
OUTPUT_PDF = REPOSITORY / "tex_source_elsevier" / "figure_8.pdf"
OUTPUT_PNG = HERE / "figure_8.png"

REDUCED_FILE = (
    VTK_ROOT
    / "reduced_xy_shear"
    / "EXAMPLE_1"
    / "vtu"
    / "EXAMPLE_1_0_1.vtu"
)

# Deformation scale for the (small-strain) shear displacements and the
# diverging colormap (signed UR_z) shared by the render and the colorbar.
WARP_FACTOR = 0.0
CMAP = "RdBu_r"

# Camera: Y up, eye in the -X/+Y/-Z octant so that +X runs to the lower left and
# +Z to the right, both into the screen. Equal magnitudes with the 35.26 deg
# elevation give the true isometric view of Figure 7 (re-oriented to Y up).
VIEW_UP = (0.0, 1.0, 0.0)
CAMERA_POSITION = (-2.0, 2.0, -2.0)

# Black mesh (element) edges drawn on every constituent.
EDGE_COLOR = (0.0, 0.0, 0.0)
EDGE_WIDTH = 0.8

# S4 shell thickness (Table 1); the flat QUAD cells are extruded to a slab of
# this thickness, rendered just outside the continuum surface.
SHELL_THICKNESS = 0.02

# elsarticle 5p/twocolumn metrics; author at the full \columnwidth so a label
# set at FS10 renders as a true 10 pt once included at ``width=\columnwidth``.
COLUMN_WIDTH_IN = 252.0 / 72.27
FS10 = 10.0 * 72.0 / 72.27
FS8 = 8.0 * 72.0 / 72.27


def read_grid(path):
    """Read an unstructured VTK grid."""
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = vtk.vtkUnstructuredGrid()
    grid.ShallowCopy(reader.GetOutput())
    return grid


def extract_cell_type(grid, cell_type):
    """Extract all cells of one VTK cell type, keeping point data."""
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
    """Minute radial offset about the origin to avoid coincident surfaces."""
    scaled = vtk.vtkUnstructuredGrid()
    scaled.DeepCopy(grid)
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(scaled.GetNumberOfPoints())
    for point_id in range(scaled.GetNumberOfPoints()):
        point = np.asarray(scaled.GetPoint(point_id)) * factor
        points.SetPoint(point_id, point)
    scaled.SetPoints(points)
    return scaled


def prepare_grid(path, warp_factor):
    """Warp the grid by ``U`` and attach the z-rotation (UR_z) scalar."""
    grid = read_grid(path)
    points = vtk_to_numpy(grid.GetPoints().GetData())
    displacement = vtk_to_numpy(grid.GetPointData().GetArray("U"))
    rotation = vtk_to_numpy(grid.GetPointData().GetArray("UR"))

    warped = points + warp_factor * displacement
    grid.GetPoints().SetData(numpy_to_vtk(np.ascontiguousarray(warped), deep=True))

    magnitude = np.linalg.norm(rotation, axis=1)
    scalar = numpy_to_vtk(np.ascontiguousarray(magnitude), deep=True)
    scalar.SetName("UR_mag")
    grid.GetPointData().AddArray(scalar)
    grid.GetPointData().SetActiveScalars("UR_mag")
    return grid, float(magnitude.max())


def build_lut(cmap, vmax, n=256):
    """Sample a Matplotlib colormap into a VTK lookup table over [0, vmax]."""
    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(n)
    lut.SetRange(0.0, vmax)
    for i in range(n):
        r, g, b, _ = cmap(i / (n - 1))
        lut.SetTableValue(i, r, g, b, 1.0)
    lut.Build()
    return lut


def color_by_scalar(mapper, lut, vmax):
    """Color a mapper by the UR_mag point scalar through the shared LUT."""
    mapper.SetLookupTable(lut)
    mapper.SetScalarRange(0.0, vmax)
    mapper.SetScalarModeToUsePointData()
    mapper.SelectColorArray("UR_mag")
    mapper.SetColorModeToMapScalars()
    mapper.ScalarVisibilityOn()


def surface_actor(grid, lut, vmax, specular, ambient, diffuse):
    """Smooth-shaded surface actor colored by the rotation magnitude."""
    surface = vtk.vtkDataSetSurfaceFilter()
    surface.SetInputData(grid)

    normals = vtk.vtkPolyDataNormals()
    normals.SetInputConnection(surface.GetOutputPort())
    normals.SetFeatureAngle(60.0)
    normals.SplittingOff()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(normals.GetOutputPort())
    color_by_scalar(mapper, lut, vmax)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    prop = actor.GetProperty()
    prop.SetInterpolationToPhong()
    prop.SetAmbient(ambient)
    prop.SetDiffuse(diffuse)
    prop.SetSpecular(specular)
    prop.SetSpecularPower(20.0)
    prop.EdgeVisibilityOn()
    prop.SetEdgeColor(*EDGE_COLOR)
    prop.SetLineWidth(EDGE_WIDTH)
    return actor


def extrude_shell_surface(grid, thickness, center):
    """Extrude the flat shell quads outward into a slab of the given thickness.

    Point normals are split at the cube edges and oriented away from the cell
    center, so each face extrudes along its own outward normal and the point
    scalar (UR_z) is carried onto the slab.
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


def shell_slab_actor(grid, lut, vmax, thickness, center,
                     specular, ambient, diffuse):
    """Smooth-shaded shell slab colored by the rotation (no element edges;
    the quad mesh is drawn separately by ``shell_edge_actor`` so the capping
    triangulation stays invisible)."""
    slab = extrude_shell_surface(grid, thickness, center)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(slab)
    color_by_scalar(mapper, lut, vmax)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    prop = actor.GetProperty()
    prop.SetInterpolationToPhong()
    prop.SetAmbient(ambient)
    prop.SetDiffuse(diffuse)
    prop.SetSpecular(specular)
    prop.SetSpecularPower(20.0)
    prop.EdgeVisibilityOff()
    return actor


def shell_edge_actor(grid, thickness, center, offset=0.0015):
    """Clean quad element edges on the outer face of the shell slab."""
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
    shifted = points + (thickness + offset) * point_normals
    poly.GetPoints().SetData(numpy_to_vtk(np.ascontiguousarray(shifted), deep=True))

    edges = vtk.vtkExtractEdges()
    edges.SetInputData(poly)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(edges.GetOutputPort())
    mapper.ScalarVisibilityOff()

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(*EDGE_COLOR)
    actor.GetProperty().SetLineWidth(EDGE_WIDTH)
    actor.GetProperty().LightingOff()
    return actor


def get_beam_radius(grid):
    """Return the beam section radius stored in the VTU, or a default."""
    radius = 0.04
    radius_array = grid.GetCellData().GetArray("RADIUS")
    if radius_array is not None and radius_array.GetNumberOfTuples() > 0:
        radius = float(radius_array.GetRange()[1])
    return radius


def beam_actor(grid, lut, vmax, radius):
    """Render beam centerlines as smooth tubes colored by the rotation."""
    geometry = vtk.vtkGeometryFilter()
    geometry.SetInputData(grid)

    tube = vtk.vtkTubeFilter()
    tube.SetInputConnection(geometry.GetOutputPort())
    tube.SetRadius(radius)
    tube.SetNumberOfSides(24)
    tube.CappingOn()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(tube.GetOutputPort())
    color_by_scalar(mapper, lut, vmax)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    prop = actor.GetProperty()
    prop.SetInterpolationToPhong()
    prop.SetAmbient(0.50)
    prop.SetDiffuse(0.60)
    prop.SetSpecular(0.0)
    prop.SetSpecularPower(20.0)
    prop.EdgeVisibilityOff()
    return actor


def beam_ring_actor(grid, radius, n_ring=24):
    """Black rings around the tubes at each finite element node (mesh edges)."""
    points = vtk.vtkPoints()
    lines = vtk.vtkCellArray()
    ring_radius = radius * 1.015
    for cell_id in range(grid.GetNumberOfCells()):
        cell = grid.GetCell(cell_id)
        if cell.GetCellType() != vtk.VTK_LINE:
            continue
        p0 = np.asarray(grid.GetPoint(cell.GetPointId(0)))
        p1 = np.asarray(grid.GetPoint(cell.GetPointId(1)))
        axis = p1 - p0
        length = np.linalg.norm(axis)
        if length < 1e-12:
            continue
        axis = axis / length
        helper = np.array([1.0, 0.0, 0.0]) if abs(axis[0]) < 0.9 \
            else np.array([0.0, 1.0, 0.0])
        u = np.cross(axis, helper)
        u /= np.linalg.norm(u)
        v = np.cross(axis, u)
        for center in (p0, p1):
            start = points.GetNumberOfPoints()
            for k in range(n_ring):
                angle = 2.0 * np.pi * k / n_ring
                q = center + ring_radius * (np.cos(angle) * u + np.sin(angle) * v)
                points.InsertNextPoint(*q)
            polyline = vtk.vtkPolyLine()
            polyline.GetPointIds().SetNumberOfIds(n_ring + 1)
            for k in range(n_ring):
                polyline.GetPointIds().SetId(k, start + k)
            polyline.GetPointIds().SetId(n_ring, start)
            lines.InsertNextCell(polyline)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetLines(lines)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polydata)
    mapper.ScalarVisibilityOff()

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(*EDGE_COLOR)
    actor.GetProperty().SetLineWidth(EDGE_WIDTH)
    actor.GetProperty().LightingOff()
    return actor


def render_cell(path, image_size=(1500, 1400)):
    """Render the warped reduced cell from a fixed isometric camera."""
    grid, vmax = prepare_grid(path, WARP_FACTOR)
    cmap = plt.get_cmap(CMAP)
    lut = build_lut(cmap, vmax)

    continuum = extract_cell_type(grid, vtk.VTK_HEXAHEDRON)
    shells = radially_scale(extract_cell_type(grid, vtk.VTK_QUAD), 1.004)
    beams = extract_cell_type(grid, vtk.VTK_LINE)

    renderer = vtk.vtkRenderer()
    renderer.SetBackground(1.0, 1.0, 1.0)
    renderer.SetBackgroundAlpha(0.0)

    light_kit = vtk.vtkLightKit()
    light_kit.AddLightsToRenderer(renderer)

    bounds = continuum.GetBounds()
    center = (0.5 * (bounds[0] + bounds[1]),
              0.5 * (bounds[2] + bounds[3]),
              0.5 * (bounds[4] + bounds[5]))
    radius = get_beam_radius(beams)
    renderer.AddActor(surface_actor(continuum, lut, vmax,
                                    specular=0.0, ambient=0.55, diffuse=0.55))
    renderer.AddActor(shell_slab_actor(shells, lut, vmax, SHELL_THICKNESS,
                                       center, specular=0.0, ambient=0.55,
                                       diffuse=0.55))
    renderer.AddActor(shell_edge_actor(shells, SHELL_THICKNESS, center))
    renderer.AddActor(beam_actor(beams, lut, vmax, radius))
    renderer.AddActor(beam_ring_actor(beams, radius))

    camera = vtk.vtkCamera()
    camera.SetPosition(*CAMERA_POSITION)
    camera.SetFocalPoint(0.0, 0.0, 0.0)
    camera.SetViewUp(*VIEW_UP)
    camera.ParallelProjectionOn()
    renderer.SetActiveCamera(camera)
    renderer.ResetCamera()
    camera.Zoom(0.92)

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
    return rgba, vmax


def crop_to_content(rgba, pad=12, thr=8):
    """Trim transparent margins around the rendered cell."""
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
    """Assemble and export the single-column figure with a horizontal colorbar.

    Authored at 0.85\\columnwidth with a fixed figure width (no tight bbox), so
    the embedded PDF keeps that natural width and the text renders at a true
    10 pt when included at ``width=0.85\\columnwidth``.
    """
    configure_matplotlib()
    image, vmax = render_cell(REDUCED_FILE)
    image = crop_to_content(image)

    width_in = 0.85 * COLUMN_WIDTH_IN
    aspect = image.shape[1] / image.shape[0]

    side = 0.015 * width_in
    img_w = width_in - 2.0 * side
    img_h = img_w / aspect

    cbar_w = 0.78 * width_in
    cbar_bar_h = 0.11
    cbar_gap = 0.07          # image bottom to colorbar bar
    cbar_labels_h = 0.50     # tick labels + axis label below the bar
    top_pad = 0.02
    bottom_pad = 0.06

    height_in = (top_pad + img_h + cbar_gap + cbar_bar_h
                 + cbar_labels_h + bottom_pad)

    figure = plt.figure(figsize=(width_in, height_in))

    image_axis = figure.add_axes([
        side / width_in,
        (cbar_labels_h + cbar_bar_h + cbar_gap + bottom_pad) / height_in,
        img_w / width_in,
        img_h / height_in,
    ])
    image_axis.imshow(image)
    image_axis.set_axis_off()

    cbar_axis = figure.add_axes([
        (width_in - cbar_w) / 2.0 / width_in,
        (bottom_pad + cbar_labels_h) / height_in,
        cbar_w / width_in,
        cbar_bar_h / height_in,
    ])

    norm = Normalize(vmin=0.0, vmax=vmax)
    mappable = cm.ScalarMappable(norm=norm, cmap=plt.get_cmap(CMAP))
    mappable.set_array([])
    colorbar = figure.colorbar(mappable, cax=cbar_axis, orientation="horizontal")
    colorbar.set_ticks(np.linspace(0.0, vmax, 5))
    colorbar.ax.xaxis.set_major_formatter(FormatStrFormatter("%.3f"))
    colorbar.ax.tick_params(labelsize=FS10, length=3.0, width=0.6, pad=2.0)
    colorbar.outline.set_linewidth(0.6)
    colorbar.set_label(r"$\|\boldsymbol{\theta}\|$", fontsize=FS10, labelpad=2.0)

    # Crop residual whitespace at the bottom while keeping the full authored
    # width, so the embedded PDF stays \columnwidth and the text true 10 pt.
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    tight = figure.get_tightbbox(renderer)
    bottom = max(tight.y0 - bottom_pad, 0.0)
    crop = Bbox([[0.0, bottom], [width_in, height_in]])

    figure.savefig(OUTPUT_PDF, dpi=600, transparent=True, bbox_inches=crop)
    figure.savefig(OUTPUT_PNG, dpi=600, bbox_inches=crop)
    plt.close(figure)

    print("Wrote {} (|UR|max = {:.6f}, warp = {:g})".format(
        OUTPUT_PDF, vmax, WARP_FACTOR))
    print("Wrote {}".format(OUTPUT_PNG))


if __name__ == "__main__":
    create_figure()
