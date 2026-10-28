"""Plot the proposed constant-rate cyclic macroscopic strain history."""

import os

import matplotlib.pyplot as plt
import numpy as np


here = os.path.dirname(os.path.abspath(__file__))
result_directory = os.path.join(here, "results")
if not os.path.isdir(result_directory):
    os.makedirs(result_directory)

peak_strains = np.asarray([
    0.0, 0.03, -0.03, 0.03, -0.03, 0.03, 0.0])
strain_jumps = np.abs(np.diff(peak_strains))
peak_times = np.r_[0.0, np.cumsum(strain_jumps)/np.sum(strain_jumps)]

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
    "mathtext.fontset": "stix",
    "font.size": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.linewidth": 0.7,
    "lines.linewidth": 1.2,
})

figure, axis = plt.subplots(figsize=(3.35, 2.35))
axis.plot(peak_times, peak_strains, color="black", marker="o",
          markersize=3.5)
axis.axhline(0.0, color="0.75", linewidth=0.6)
axis.set_xlabel(r"Pseudo-time $\tau$")
axis.set_ylabel(r"Macroscopic strain $\overline{\epsilon}_{11}$")
axis.set_xlim(0.0, 1.0)
axis.tick_params(direction="in", top=True, right=True)
figure.tight_layout()

for extension, dpi in (("pdf", None), ("png", 300)):
    figure.savefig(
        os.path.join(result_directory,
                     "cyclic_load_path_preview." + extension),
        dpi=dpi, bbox_inches="tight")

print("Peak times: {}".format(peak_times))
print("Segment slopes: {}".format(
    np.diff(peak_strains)/np.diff(peak_times)))
