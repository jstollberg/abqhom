"""Plot the initial and evolved effective yield surfaces."""

import csv
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


here = os.path.dirname(os.path.abspath(__file__))
yield_surface_directory = os.path.join(here, "yield_surface")
surface_directory = os.path.join(yield_surface_directory, "surface_3d")

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
    "mathtext.fontset": "stix",
    "font.size": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.linewidth": 0.7,
})


def read_surface(state):
    """Return normal stresses, base stress, and triangular faces."""
    directory = os.path.join(surface_directory, state)
    with open(os.path.join(directory, "yield_surface.csv"), newline="") \
            as csv_file:
        rows = list(csv.DictReader(csv_file))
    stress = np.asarray([
        [float(row["sigma-11"]), float(row["sigma-22"]),
         float(row["sigma-33"])] for row in rows])
    base = np.asarray([
        [float(row["base-sigma-11"]), float(row["base-sigma-22"]),
         float(row["base-sigma-33"])] for row in rows])
    faces = np.loadtxt(os.path.join(directory, "faces.csv"), delimiter=",",
                       skiprows=1, dtype=int)
    return stress, base[0], faces


def add_surface(axis, points, faces, color, alpha=0.55):
    """Add one triangulated surface to a three-dimensional axis."""
    collection = Poly3DCollection(
        points[faces], facecolor=color, alpha=alpha,
        edgecolor="none", linewidth=0.0)
    collection.set_rasterized(True)
    axis.add_collection3d(collection)


def convex_hull(points):
    """Return the counterclockwise convex hull of two-dimensional points."""
    points = np.unique(points, axis=0)
    points = points[np.lexsort((points[:, 1], points[:, 0]))]

    def cross(origin, first, second):
        first = first - origin
        second = second - origin
        return first[0]*second[1] - first[1]*second[0]

    lower = []
    for point in points:
        while len(lower) >= 2 and cross(
                lower[-2], lower[-1], point) <= 0.0:
            lower.pop()
        lower.append(point)
    upper = []
    for point in points[::-1]:
        while len(upper) >= 2 and cross(
                upper[-2], upper[-1], point) <= 0.0:
            upper.pop()
        upper.append(point)
    return np.asarray(lower[:-1] + upper[:-1])


def add_projection(axis, initial, evolved, first, second):
    """Plot the projected outlines of both yield surfaces."""
    initial_hull = convex_hull(initial[:, (first, second)])
    evolved_hull = convex_hull(evolved[:, (first, second)])
    for hull, color, linestyle in (
            (initial_hull, "#0072B2", "--"),
            (evolved_hull, "#D55E00", "-")):
        closed = np.vstack((hull, hull[0]))
        axis.fill(hull[:, 0], hull[:, 1], color=color, alpha=0.12)
        axis.plot(closed[:, 0], closed[:, 1], color=color,
                  linestyle=linestyle, linewidth=1.2)
    axis.scatter(base_transformed[first], base_transformed[second],
                 marker="x", s=20, linewidth=1.0, color="black", zorder=3)
    axis.axhline(0.0, color="0.8", linewidth=0.5, zorder=0)
    axis.axvline(0.0, color="0.8", linewidth=0.5, zorder=0)
    axis.set_xlabel(labels[first])
    axis.set_ylabel(labels[second])
    axis.set_aspect("equal", adjustable="box")
    axis.tick_params(direction="in", top=True, right=True)


initial, _, initial_faces = read_surface("initial")
evolved, base, evolved_faces = read_surface("final")
transformation = np.asarray([
    [1.0/np.sqrt(3.0), 1.0/np.sqrt(3.0), 1.0/np.sqrt(3.0)],
    [2.0/np.sqrt(6.0), -1.0/np.sqrt(6.0), -1.0/np.sqrt(6.0)],
    [0.0, 1.0/np.sqrt(2.0), -1.0/np.sqrt(2.0)],
])
initial_transformed = initial.dot(transformation.T)
evolved_transformed = evolved.dot(transformation.T)
base_transformed = transformation.dot(base)
limits = np.asarray([
    np.minimum(initial_transformed.min(axis=0),
               evolved_transformed.min(axis=0)),
    np.maximum(initial_transformed.max(axis=0),
               evolved_transformed.max(axis=0)),
])
margin = 0.06*(limits[1] - limits[0])
limits[0] -= margin
limits[1] += margin

labels = (
    r"$\overline{\sigma}_{\mathrm{h}}"
    r"$",
    r"$\overline{\sigma}_{\mathrm{d}1}"
    r"$",
    r"$\overline{\sigma}_{\mathrm{d}2}"
    r"$",
)

figure = plt.figure(figsize=(6.8, 5.2))
axis_3d = figure.add_subplot(2, 2, 1, projection="3d")
projection_axes = [figure.add_subplot(2, 2, index)
                   for index in (2, 3, 4)]

add_surface(axis_3d, initial_transformed, initial_faces,
            "#0072B2", alpha=0.35)
add_surface(axis_3d, evolved_transformed, evolved_faces,
            "#D55E00", alpha=0.55)
axis_3d.scatter(*base_transformed, marker="x", s=20, linewidth=1.0,
                color="black", depthshade=False)
axis_3d.set_xlim(*limits[:, 0])
axis_3d.set_ylim(*limits[:, 1])
axis_3d.set_zlim(*limits[:, 2])
axis_3d.set_box_aspect(limits[1] - limits[0])
axis_3d.set_xlabel(r"$\overline{\sigma}_{\mathrm{h}}$", labelpad=-1)
axis_3d.set_ylabel(r"$\overline{\sigma}_{\mathrm{d}1}$", labelpad=-1)
axis_3d.set_zlabel(r"$\overline{\sigma}_{\mathrm{d}2}$", labelpad=-2)
axis_3d.view_init(elev=20.0, azim=-55.0)
axis_3d.set_proj_type("ortho")
axis_3d.tick_params(pad=-2)

for axis, coordinates, title in zip(
        projection_axes, ((0, 1), (0, 2), (1, 2)),
        ("Hydrostatic--deviatoric projection",
         "Hydrostatic--deviatoric projection",
         "Deviatoric projection")):
    add_projection(axis, initial_transformed, evolved_transformed,
                   *coordinates)
    axis.set_xlim(*limits[:, coordinates[0]])
    axis.set_ylim(*limits[:, coordinates[1]])
    axis.set_title(title, pad=3.0)

for panel, axis in zip(("(a)", "(b)", "(c)", "(d)"),
                       (axis_3d, *projection_axes)):
    if axis is axis_3d:
        axis.text2D(0.0, 1.0, panel, transform=axis.transAxes,
                    ha="left", va="top", fontweight="bold")
    else:
        axis.text(0.0, 1.02, panel, transform=axis.transAxes,
                  ha="left", va="bottom", fontweight="bold")

handles = (
    Patch(facecolor="#0072B2", alpha=0.35, label="Initial surface"),
    Patch(facecolor="#D55E00", alpha=0.55, label="Evolved surface"),
    Line2D([], [], color="black", marker="x", linestyle="none",
           label="Probing state"),
)
figure.legend(handles=handles, frameon=False, loc="upper center",
              bbox_to_anchor=(0.5, 1.0), ncol=3)
figure.subplots_adjust(left=0.09, right=0.98, bottom=0.1, top=0.91,
                      wspace=0.32, hspace=0.38)

output = os.path.join(surface_directory, "yield_surface_3d")
figure.savefig(output + ".pdf", bbox_inches="tight")
figure.savefig(output + ".png", dpi=300, bbox_inches="tight")

# Standalone view of the initial surface along the hydrostatic axis.
initial_hull = convex_hull(initial_transformed[:, (1, 2)])
closed_hull = np.vstack((initial_hull, initial_hull[0]))
hydrostatic_figure, hydrostatic_axis = plt.subplots(figsize=(3.4, 3.2))
hydrostatic_axis.fill(initial_hull[:, 0], initial_hull[:, 1],
                      color="#0072B2", alpha=0.18)
hydrostatic_axis.plot(closed_hull[:, 0], closed_hull[:, 1],
                      color="#0072B2", linewidth=1.4)
hydrostatic_axis.axhline(0.0, color="0.8", linewidth=0.5, zorder=0)
hydrostatic_axis.axvline(0.0, color="0.8", linewidth=0.5, zorder=0)
hydrostatic_axis.set_xlabel(labels[1])
hydrostatic_axis.set_ylabel(labels[2])
hydrostatic_axis.set_aspect("equal", adjustable="box")
hydrostatic_axis.tick_params(direction="in", top=True, right=True)
hydrostatic_axis.set_title(r"Initial surface viewed along $\overline{\sigma}_{\mathrm{h}}$",
                           pad=5.0)
hydrostatic_figure.tight_layout()

hydrostatic_output = os.path.join(
    surface_directory, "initial_surface_hydrostatic_axis")
hydrostatic_figure.savefig(hydrostatic_output + ".pdf",
                           bbox_inches="tight")
hydrostatic_figure.savefig(hydrostatic_output + ".png", dpi=300,
                           bbox_inches="tight")

# Conventional principal-stress-space representation.  The hydrostatic axis
# is the diagonal sigma_11 = sigma_22 = sigma_33.
principal_figure = plt.figure(figsize=(4.2, 3.8))
principal_axis = principal_figure.add_subplot(111, projection="3d")
add_surface(principal_axis, initial, initial_faces, "#0072B2", alpha=0.58)

stress_limit = 1.08*np.max(np.abs(initial))
hydrostatic_coordinate = np.asarray([-stress_limit, stress_limit])
principal_axis.plot(hydrostatic_coordinate, hydrostatic_coordinate,
                    hydrostatic_coordinate, color="black", linestyle="--",
                    linewidth=1.2,
                    label=r"Hydrostatic axis: $\overline{\sigma}_{11}"
                          r"=\overline{\sigma}_{22}"
                          r"=\overline{\sigma}_{33}$")
principal_axis.scatter(0.0, 0.0, 0.0, color="black", s=8,
                       depthshade=False)
principal_axis.set_xlim(-stress_limit, stress_limit)
principal_axis.set_ylim(-stress_limit, stress_limit)
principal_axis.set_zlim(-stress_limit, stress_limit)
principal_axis.set_box_aspect((1.0, 1.0, 1.0))
principal_axis.set_xlabel(r"$\overline{\sigma}_{11}$", labelpad=0)
principal_axis.set_ylabel(r"$\overline{\sigma}_{22}$", labelpad=0)
principal_axis.set_zlabel(r"$\overline{\sigma}_{33}$", labelpad=0)
principal_axis.view_init(elev=22.0, azim=-48.0)
principal_axis.set_proj_type("ortho")
principal_axis.tick_params(pad=-1)
principal_axis.legend(frameon=False, loc="upper center",
                      bbox_to_anchor=(0.5, 1.04), fontsize=7)
principal_figure.subplots_adjust(left=0.0, right=0.94, bottom=0.02, top=0.9)

principal_output = os.path.join(
    surface_directory, "initial_surface_principal_stress_axes")
principal_figure.savefig(principal_output + ".pdf", bbox_inches="tight")
principal_figure.savefig(principal_output + ".png", dpi=300,
                         bbox_inches="tight")

# View directly along the hydrostatic direction.  The projections of the
# three normal-stress axes have equal length and enclose angles of 120 deg.
isometric_hull = convex_hull(initial_transformed[:, (1, 2)])
isometric_closed = np.vstack((isometric_hull, isometric_hull[0]))
isometric_figure, isometric_axis = plt.subplots(figsize=(3.6, 3.4))
isometric_axis.fill(isometric_hull[:, 0], isometric_hull[:, 1],
                    color="#0072B2", alpha=0.2, zorder=1)
isometric_axis.plot(isometric_closed[:, 0], isometric_closed[:, 1],
                    color="#0072B2", linewidth=1.4, zorder=2)

projected_axes = np.asarray([
    [2.0/np.sqrt(6.0), 0.0],
    [-1.0/np.sqrt(6.0), 1.0/np.sqrt(2.0)],
    [-1.0/np.sqrt(6.0), -1.0/np.sqrt(2.0)],
])
axis_radius = 1.25*np.max(np.linalg.norm(isometric_hull, axis=1))
projected_axes /= np.linalg.norm(projected_axes, axis=1)[:, None]
axis_labels = (r"$\overline{\sigma}_{11}$",
               r"$\overline{\sigma}_{22}$",
               r"$\overline{\sigma}_{33}$")
for direction, label in zip(projected_axes, axis_labels):
    endpoint = axis_radius*direction
    isometric_axis.annotate(
        "", xy=endpoint, xytext=(0.0, 0.0),
        arrowprops={"arrowstyle": "->", "color": "0.25",
                    "linewidth": 0.8}, zorder=0)
    label_position = 1.09*endpoint
    isometric_axis.text(label_position[0], label_position[1], label,
                        ha="center", va="center")

# A circled dot denotes an axis directed out of the page and toward the
# observer.
symbol_position = np.asarray([-0.94*axis_radius, 0.94*axis_radius])
isometric_axis.scatter(*symbol_position, s=42, facecolors="none",
                       edgecolors="black", linewidths=0.8, zorder=4)
isometric_axis.scatter(*symbol_position, s=5, color="black", zorder=5)
isometric_axis.text(symbol_position[0] + 0.08*axis_radius,
                    symbol_position[1],
                    r"$\overline{\sigma}_{\mathrm{h}}$",
                    ha="left", va="center")

plot_limit = 1.43*axis_radius
isometric_axis.set_xlim(-plot_limit, plot_limit)
isometric_axis.set_ylim(-plot_limit, plot_limit)
isometric_axis.set_aspect("equal", adjustable="box")
isometric_axis.set_axis_off()
isometric_axis.set_title("Initial surface viewed along the hydrostatic axis",
                         pad=3.0)
isometric_figure.tight_layout()

isometric_output = os.path.join(
    surface_directory, "initial_surface_hydrostatic_isometric")
isometric_figure.savefig(isometric_output + ".pdf", bbox_inches="tight")
isometric_figure.savefig(isometric_output + ".png", dpi=300,
                         bbox_inches="tight")
