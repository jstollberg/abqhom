import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from abqhom.utils import find_all_files, read_csv_file

base_path = os.path.join("C:/Users/jonat/Documents/abqhom/examples",
                         "fiber_reinforced_2d/results")
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
for d, p, r, v in zip(files_drlt, files_periodic, files_reuss, files_voigt):
    drlt.append(read_csv_file(d))
    periodic.append(read_csv_file(p))
    reuss.append(read_csv_file(r))
    voigt.append(read_csv_file(v))

    file_name = os.path.basename(d)
    volfrac, _ = os.path.splitext(file_name)
    volfracs.append(float(volfrac))

drlt = np.array(drlt)
periodic = np.array(periodic)
reuss = np.array(reuss)
voigt = np.array(voigt)
volfracs = np.array(volfracs)

# Reorder eigenvalues
eig3 = np.vstack((drlt[:, 0], periodic[:, 0], reuss[:, 0], voigt[:, 0]))
eig2 = np.vstack((drlt[:, 1], periodic[:, 1], reuss[:, 1], voigt[:, 1]))
eig1 = np.vstack((drlt[:, 2], periodic[:, 2], reuss[:, 2], voigt[:, 2]))

# Set font size
matplotlib.rcParams.update({"font.size": 14})

# Plot first eigenvalue
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(36/2.54, 15/2.54),
                                    dpi=300)
ax1.plot(volfracs, eig1[2], "--", color="black", label="Reuss lower bound",
         linewidth=3)
ax1.plot(volfracs, eig1[3], ":", color="black", label="Voigt upper bound",
         linewidth=3)
ax1.plot(volfracs, eig1[0], marker=".", color="blue", label="Numerical (ADBC)",
         linewidth=3, markersize=12)
ax1.plot(volfracs, eig1[1], marker=".", color="red", label="Numerical (PBC)",
         linewidth=3, markersize=12)
ax1.set(xlabel="Fiber volume fraction", ylabel="First eigenvalue")
ax1.legend()

# Plot second eigenvalue
ax2.plot(volfracs, eig2[2], "--", color="black", label="Reuss lower bound",
         linewidth=3)
ax2.plot(volfracs, eig2[3], ":", color="black", label="Voigt upper bound",
         linewidth=3)
ax2.plot(volfracs, eig2[0], marker=".", color="blue", label="Numerical (ADBC)",
         linewidth=3, markersize=12)
ax2.plot(volfracs, eig2[1], marker=".", color="red", label="Numerical (PBC)",
         linewidth=3, markersize=12)
ax2.set(xlabel="Fiber volume fraction", ylabel="Second eigenvalue")

# Plot third eigenvalue
ax3.plot(volfracs, eig3[2], "--", color="black", label="Reuss lower bound",
         linewidth=3)
ax3.plot(volfracs, eig3[3], ":", color="black", label="Voigt upper bound",
         linewidth=3)
ax3.plot(volfracs, eig3[0], marker=".", color="blue", label="Numerical (ADBC)",
         linewidth=3, markersize=12)
ax3.plot(volfracs, eig3[1], marker=".", color="red", label="Numerical (PBC)",
         linewidth=3, markersize=12)
ax3.set(xlabel="Fiber volume fraction", ylabel="Third eigenvalue")

fig.tight_layout()
fig.savefig("eigenvalues.png")
plt.show()
