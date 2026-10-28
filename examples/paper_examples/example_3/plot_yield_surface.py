"""Plot the initial and evolved effective yield-surface sections."""

import csv
import os

import matplotlib.pyplot as plt
import numpy as np


here = os.path.dirname(os.path.abspath(__file__))
yield_surface_directory = os.path.join(here, "yield_surface")
initial_directory = os.path.join(yield_surface_directory, "results")
final_directory = os.path.join(yield_surface_directory, "final_results")

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
    "mathtext.fontset": "stix",
    "font.size": 8,
    "axes.labelsize": 8,
    "legend.fontsize": 7,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.linewidth": 0.7,
    "lines.linewidth": 1.0,
})


def read_results(directory, name):
    """Read one yield-surface section, including its phase labels."""
    with open(os.path.join(directory, name + ".csv"), newline="") \
            as csv_file:
        rows = list(csv.DictReader(csv_file))
    numeric_names = [name for name in rows[0]
                     if name != "controlling-phase"]
    numeric = {name: np.asarray([float(row[name]) for row in rows])
               for name in numeric_names}
    numeric["controlling-phase"] = np.asarray(
        [row["controlling-phase"] for row in rows])
    return numeric


sections = (
    ("sigma11_sigma22", "sigma-11", "sigma-22",
     r"$\overline{\sigma}_{11}$", r"$\overline{\sigma}_{22}$"),
    ("sigma11_sigma33", "sigma-11", "sigma-33",
     r"$\overline{\sigma}_{11}$", r"$\overline{\sigma}_{33}$"),
)
phase_style = {
    "SHELL": ("o", "#009E73", "Shell faces yield first"),
    "BEAM": ("s", "#CC79A7", "Beam core yields first"),
    "BOTH": ("D", "#0072B2", "Simultaneous onset"),
}

figure, axes = plt.subplots(1, 2, figsize=(4.8, 2.4))
used_labels = set()
for panel, axis, section in zip(("(a)", "(b)"), axes, sections):
    name, x_name, y_name, x_label, y_label = section
    initial = read_results(initial_directory, name)
    initial_x = initial[x_name]
    initial_y = initial[y_name]
    axis.plot(np.r_[initial_x, initial_x[0]],
              np.r_[initial_y, initial_y[0]],
              color="black",
              label="Initial surface" if panel == "(a)" else None,
              zorder=1)
    x = initial_x
    y = initial_y
    for phase, (marker, color, label) in phase_style.items():
        selected = initial["controlling-phase"] == phase
        if not np.any(selected):
            continue
        display_label = label if label not in used_labels else None
        axis.scatter(x[selected], y[selected], marker=marker, s=15,
                     facecolor="white", edgecolor=color, linewidth=0.9,
                     label=display_label, zorder=2)
        used_labels.add(label)
    axis.axhline(0.0, color="0.75", linewidth=0.6)
    axis.axvline(0.0, color="0.75", linewidth=0.6)
    axis.set_xlabel(x_label)
    axis.set_ylabel(y_label)
    axis.set_aspect("equal", adjustable="box")
    axis.tick_params(direction="in", top=True, right=True)
    axis.text(0.0, 1.04, panel, transform=axis.transAxes,
              ha="left", va="bottom", fontweight="bold")

handles = []
labels = []
for axis in axes:
    axis_handles, axis_labels = axis.get_legend_handles_labels()
    handles.extend(axis_handles)
    labels.extend(axis_labels)
if handles:
    figure.legend(handles, labels, frameon=False, loc="upper center",
                  bbox_to_anchor=(0.5, 1.08), ncol=len(handles))

figure.tight_layout(w_pad=1.5)
for extension, dpi in (("pdf", None), ("png", 300)):
    figure.savefig(os.path.join(
        initial_directory, "initial_yield_surface." + extension),
        dpi=dpi, bbox_inches="tight")


# Comparison after monotonic loading and strictly elastic unloading to the
# common radial-probing state.
comparison, axes = plt.subplots(1, 2, figsize=(4.8, 2.4))
for panel, axis, section in zip(("(a)", "(b)"), axes, sections):
    name, x_name, y_name, x_label, y_label = section
    initial = read_results(initial_directory, name)
    evolved = read_results(final_directory, name)

    initial_x = initial[x_name]
    initial_y = initial[y_name]
    evolved_x = evolved[x_name]
    evolved_y = evolved[y_name]
    axis.plot(np.r_[initial_x, initial_x[0]],
              np.r_[initial_y, initial_y[0]],
              color="0.55", linestyle="--",
              label="Initial surface" if panel == "(a)" else None,
              zorder=1)
    axis.plot(np.r_[evolved_x, evolved_x[0]],
              np.r_[evolved_y, evolved_y[0]],
              color="black",
              label="Evolved surface" if panel == "(a)" else None,
              zorder=2)
    axis.scatter(evolved["base-sigma-11"][0],
                 evolved["base-sigma-22" if name.endswith("22")
                         else "base-sigma-33"][0],
                 marker="x", s=24, linewidth=1.0, color="#D55E00",
                 label="Probing state" if panel == "(a)" else None,
                 zorder=3)
    axis.axhline(0.0, color="0.75", linewidth=0.6)
    axis.axvline(0.0, color="0.75", linewidth=0.6)
    axis.set_xlabel(x_label)
    axis.set_ylabel(y_label)
    axis.set_aspect("equal", adjustable="box")
    axis.tick_params(direction="in", top=True, right=True)
    axis.text(0.0, 1.04, panel, transform=axis.transAxes,
              ha="left", va="bottom", fontweight="bold")

handles, labels = axes[0].get_legend_handles_labels()
comparison.legend(handles, labels, frameon=False, loc="upper center",
                  bbox_to_anchor=(0.5, 1.08), ncol=len(handles))
comparison.tight_layout(w_pad=1.5)
for extension, dpi in (("pdf", None), ("png", 300)):
    comparison.savefig(os.path.join(
        final_directory, "yield_surface_evolution." + extension),
        dpi=dpi, bbox_inches="tight")

