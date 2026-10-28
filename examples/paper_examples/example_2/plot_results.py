"""Plot engineering deformation against effective Cauchy stress."""

import os

import matplotlib.pyplot as plt
import numpy as np


here = os.path.dirname(os.path.abspath(__file__))
result_root = os.path.join(here, "final_results")
cases = (
    ("compression", "Compression", r"Compression $1-\lambda_1$", -1.0),
    ("tension", "Tension", r"Extension $\lambda_1-1$", 1.0),
    ("shear_12", "Shear", r"Shear $\overline{H}_{12}$", 1.0),
)

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 9,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.5,
})

figure, axes = plt.subplots(1, 3, figsize=(7.2, 2.35))
for axis, (case, title, xlabel, sign) in zip(axes, cases):
    response = np.loadtxt(
        os.path.join(result_root, case, "stress_curve.csv"), delimiter=",")
    if case == "compression":
        abscissa = 1.0 - response[:, 0]
    elif case == "tension":
        abscissa = response[:, 0] - 1.0
    else:
        abscissa = response[:, 0]
    axis.plot(abscissa, sign*response[:, 2], color="#1f77b4")
    axis.set_title(title)
    axis.set_xlabel(xlabel)
    axis.set_xlim(0.0, 0.4)
    axis.set_ylim(bottom=0.0)
    axis.grid(alpha=0.2, linewidth=0.5)
    axis.tick_params(direction="in")

axes[0].set_ylabel(r"Effective Cauchy stress magnitude")
figure.tight_layout(w_pad=1.2)
output = os.path.join(result_root, "effective_response.png")
figure.savefig(output, dpi=300, bbox_inches="tight")
print(output)
