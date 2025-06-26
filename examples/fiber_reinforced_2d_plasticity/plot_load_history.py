import numpy as np
import matplotlib
import matplotlib.pyplot as plt


def load_path(time):
    strain = np.zeros(6)
    strain[0] = np.sin(np.pi*time)
    strain[1] = np.exp(-((time - 0.5)/0.2)**2) - np.exp(-(0.5/0.2)**2)
    strain[5] = np.sin(2*np.pi*time)
    strain *= 0.125
    return strain


# Compute data points
time = np.linspace(0.0, 1.0, 50)
strain = []
for t in time:
    strain.append(load_path(t))
strain = np.array(strain)

# Font settings
matplotlib.rcParams.update({"font.size": 15})
plt.rcParams['text.usetex'] = True

# Plot the strain history
fig, ax = plt.subplots(figsize=(18/2.54, 12/2.54), dpi=300)
ax.plot(time, strain[:, 0], label=r"$\varepsilon_{11}$", linewidth=2,
        color="black")
ax.plot(time, strain[:, 1], "--", label=r"$\varepsilon_{22}$", linewidth=2,
        color="black")
ax.plot(time, strain[:, 5], ":", label=r"$\varepsilon_{12}$", linewidth=2,
        color="black")

# Axis settings
ax.set(xlabel="Pseudo-time", ylabel="Strain")
ax.set_xbound(0.0, 1.0)
ax.legend()
ax.grid()

fig.tight_layout()
fig.savefig("load_history.png")
plt.show()
