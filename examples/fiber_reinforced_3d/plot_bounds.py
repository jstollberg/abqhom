import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from abqhom.utils import find_all_files, read_csv_file

base_path = os.path.join("C:/Users/jonat/Documents/abqhom/examples",
                         "fiber_reinforced_3d/results")
drlt_path = os.path.join(base_path, "DIRICHLET")
periodic_path = os.path.join(base_path, "PERIODIC")

# Find files
files_drlt = find_all_files(os.path.join(drlt_path, "numerical"),
                            extension=".csv")
files_periodic = find_all_files(os.path.join(periodic_path, "numerical"),
                                extension=".csv")
files_reuss = find_all_files(os.path.join(drlt_path, "reuss"),
                             extension=".csv")
files_voigt = find_all_files(os.path.join(drlt_path, "voigt"),
                             extension=".csv")

# Collect data
drlt = []
periodic = []
reuss = []
voigt = []
volfracs = []
for i in range(min(len(files_drlt), len(files_periodic),
                   len(files_reuss), len(files_voigt))):
    drlt.append(read_csv_file(files_drlt[i]))
    periodic.append(read_csv_file(files_periodic[i]))
    reuss.append(read_csv_file(files_reuss[i]))
    voigt.append(read_csv_file(files_voigt[i]))

    file_name = os.path.basename(files_drlt[i])
    volfrac, _ = os.path.splitext(file_name)
    volfracs.append(float(volfrac))

drlt = np.array(drlt)
periodic = np.array(periodic)
reuss = np.array(reuss)
voigt = np.array(voigt)
volfracs = np.array(volfracs)

# Extract eigenvalues
eigs = []
for i in range(6):
    eigs.append(np.vstack((drlt[:, i], periodic[:, i],
                           reuss[:, i], voigt[:, i])))
eigs = eigs[::-1]

# Set font size
matplotlib.rcParams.update({"font.size": 14})

# Plot first eigenvalue
labels = ["First", "Second", "Third", "Fourth", "Fifth", "Sixth"]
fig, axes = plt.subplots(2, 3, figsize=(36/2.54, 30/2.54), dpi=300)

for ax, eig, label in zip(axes.flatten(), eigs, labels):
    ax.plot(volfracs, eig[2], "--", color="black", label="Reuss lower bound",
            linewidth=3)
    ax.plot(volfracs, eig[3], ":", color="black", label="Voigt upper bound",
            linewidth=3)
    ax.plot(volfracs, eig[0], marker=".", color="blue",
            label="Numerical (ADBC)", linewidth=3, markersize=12)
    ax.plot(volfracs, eig[1], marker=".", color="red", label="Numerical (PBC)",
            linewidth=3, markersize=12)
    ax.set(xlabel="Fiber volume fraction", ylabel=label + " eigenvalue")

    if label == "First":
        ax.legend()

fig.tight_layout()
fig.savefig("reuss_voigt_bounds.png")
plt.show()
