"""Plot the macroscopic response and phase-wise plastic-strain history."""

import os

import matplotlib.pyplot as plt
import numpy as np


here = os.path.dirname(os.path.abspath(__file__))
result_directory = os.path.join(here, "results")

response = np.genfromtxt(
    os.path.join(result_directory, "cyclic_response.csv"),
    delimiter=",", names=True)
plasticity = np.genfromtxt(
    os.path.join(result_directory, "phase_plasticity.csv"),
    delimiter=",", names=True)

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
    "lines.linewidth": 1.2,
})

stress_data = (
    ("sigma_11", r"$\overline{\sigma}_{11}$", "black"),
    ("sigma_22", r"$\overline{\sigma}_{22}$", "#0072B2"),
    ("sigma_33", r"$\overline{\sigma}_{33}$", "#D55E00"),
)
phase_data = (
    ("max_peeq_shell", "Shell faces", "#009E73"),
    ("max_peeq_beam", "Beam core", "#CC79A7"),
)


def finish_axes(axes):
    """Apply the common axis and panel-label formatting."""
    for panel, axis in zip(("(a)", "(b)", "(c)"), axes):
        axis.tick_params(direction="in", top=True, right=True)
        axis.text(0.0, 1.04, panel, transform=axis.transAxes,
                  ha="left", va="bottom", fontweight="bold")


def save_figure(figure, name):
    """Save a figure in vector and raster formats."""
    for extension, dpi in (("pdf", None), ("png", 300)):
        figure.savefig(
            os.path.join(result_directory, name + "." + extension),
            dpi=dpi, bbox_inches="tight")


# Standalone normal stress--strain response.
stress_figure, stress_axes = plt.subplots(
    1, 3, figsize=(7.1, 2.35), sharex=True)
for axis, (field, label, color) in zip(stress_axes, stress_data):
    axis.plot(response["epsilon_11"], response[field], color=color)
    axis.axhline(0.0, color="0.75", linewidth=0.6)
    axis.axvline(0.0, color="0.75", linewidth=0.6)
    axis.set_xlabel(r"Macroscopic strain $\overline{\epsilon}_{11}$")
    axis.set_title(label, pad=5.0)
stress_axes[0].set_ylabel("Effective stress")
finish_axes(stress_axes)
stress_figure.tight_layout(w_pad=1.0)
save_figure(stress_figure, "stress_strain_response")


# Supporting histories of the prescribed loading and phase plasticity.
history_figure, history_axes = plt.subplots(1, 2, figsize=(4.8, 2.35))
history_axis, plasticity_axis = history_axes
history_axis.plot(response["time"], response["epsilon_11"], color="black")
history_axis.axhline(0.0, color="0.75", linewidth=0.6)
history_axis.set_xlabel(r"Pseudo-time $\tau$")
history_axis.set_ylabel(r"Macroscopic strain $\overline{\epsilon}_{11}$")
history_axis.set_xlim(response["time"][0], response["time"][-1])

for field, label, color in phase_data:
    plasticity_axis.plot(plasticity["time"], plasticity[field],
                         label=label, color=color)
plasticity_axis.set_xlabel(r"Pseudo-time $\tau$")
plasticity_axis.set_ylabel(r"Maximum equivalent plastic strain")
plasticity_axis.set_xlim(plasticity["time"][0], plasticity["time"][-1])
plasticity_axis.set_ylim(bottom=0.0)
plasticity_axis.legend(frameon=False, loc="upper left")

for panel, axis in zip(("(a)", "(b)"), history_axes):
    axis.tick_params(direction="in", top=True, right=True)
    axis.text(0.0, 1.04, panel, transform=axis.transAxes,
              ha="left", va="bottom", fontweight="bold")
history_figure.tight_layout(w_pad=1.5)
save_figure(history_figure, "cyclic_plasticity_response")

