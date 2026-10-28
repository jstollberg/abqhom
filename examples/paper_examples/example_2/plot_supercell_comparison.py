"""Plot the 2 x 2 x 2 supercell compression responses."""

import os

import matplotlib.pyplot as plt
import numpy as np


here = os.path.dirname(os.path.abspath(__file__))
result_root = os.path.join(here, "supercell_results")
cases = (
    ("comp_n16_mode_01", "lowest primitive-periodic mode", "#1f77b4"),
    ("comp_n16_mode_03", "lowest non-primitive mode", "#d62728"),
)

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 9,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.5,
})

figure, axis = plt.subplots(figsize=(3.45, 2.7))
for case, label, color in cases:
    response = np.loadtxt(
        os.path.join(result_root, case, "stress_curve.csv"), delimiter=",")
    compression = 1.0 - response[:, 1]
    axis.plot(compression, -response[:, 3], label=label, color=color)

axis.set_xlabel(r"Engineering compression $1-\lambda_1$")
axis.set_ylabel(r"Effective Cauchy stress magnitude")
axis.set_xlim(0.0, 0.4)
axis.set_ylim(bottom=0.0)
axis.grid(alpha=0.2, linewidth=0.5)
axis.tick_params(direction="in")
axis.legend(frameon=False, fontsize=8)
figure.tight_layout()

for extension in ("png", "pdf"):
    output = os.path.join(
        result_root, "compression_mode_comparison." + extension)
    figure.savefig(output, dpi=300, bbox_inches="tight")
    print(output)


# Direct cell-size comparison with identical discretizations and lowest-mode
# imperfections.
single_cell = np.loadtxt(
    os.path.join(here, "final_results", "compression", "stress_curve.csv"),
    delimiter=",")
supercell = np.loadtxt(
    os.path.join(result_root, "comp_n16_mode_01", "stress_curve.csv"),
    delimiter=",")

figure, axis = plt.subplots(figsize=(3.45, 2.7))
axis.plot(1.0 - single_cell[:, 0], -single_cell[:, 2],
          label=r"single cell ($1^3$)", color="#1f77b4")
axis.plot(1.0 - supercell[:, 1], -supercell[:, 3],
          label=r"supercell ($2^3$)", color="#d62728", linestyle="--")
axis.set_xlabel(r"Engineering compression $1-\lambda_1$")
axis.set_ylabel(r"Effective Cauchy stress magnitude")
axis.set_xlim(0.0, 0.4)
axis.set_ylim(bottom=0.0)
axis.grid(alpha=0.2, linewidth=0.5)
axis.tick_params(direction="in")
axis.legend(frameon=False, fontsize=8)
figure.tight_layout()

for extension in ("png", "pdf"):
    output = os.path.join(
        result_root, "compression_cell_size_comparison." + extension)
    figure.savefig(output, dpi=300, bbox_inches="tight")
    print(output)
