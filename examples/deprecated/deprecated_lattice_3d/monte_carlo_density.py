import numpy as np
import os
from abqhom.RVE import finalize_model, create_beam_structure
from abqhom.latticetools import approx_rel_density, monte_carlo_density
from abqhom.examples import simple_RVE_3d
from abqhom.utils import export_csv_file

# set working directory
base_path = "C:/Users/jonat/Documents/"
workdir = os.path.join(base_path, 
                       "Institute for Mechanics/abqhom/examples",
                       "lattice_foam_3d")
if not os.path.isdir(workdir):
    os.mkdir(workdir)
os.chdir(workdir)

# create folders to store results
monte_carlo_path = os.path.join(workdir, "monte_carlo_density_results")
if not os.path.isdir(monte_carlo_path):
    os.mkdir(monte_carlo_path)
approx_path = os.path.join(workdir, "approx_density_results")
if not os.path.isdir(approx_path):
    os.mkdir(approx_path)
    
# RVE size
x = y = z = 0.0
dx = dy = dz = 10
lc = 1.25

# monte carlo parameters
tol = 1e-3
min_points = 1000
max_points = 2000

# parameter range
aspect_ratios = np.linspace(0.01, 1.0, 50)

# create geometry in gmsh
model_name, group_map = simple_RVE_3d(dx=dx, dy=dy, dz=dz, lc=lc)
node_tags, node_coords, el_tags, conn = create_beam_structure(model_name)
    
# run integration schemes
density_approx = []
density_monte_carlo = []
for ar in aspect_ratios:
    print("aspect_radio =", ar)
    radius = (ar*lc)/2
    
    # compute relative density without considering overlapping at joints
    rho_1 = approx_rel_density(model_name, group_map, node_tags, node_coords, 
                             el_tags, conn, radius, dx, dy, dz)
    density_approx.append(rho_1)
    file = "approx_{}.csv".format(ar)
    export_csv_file(np.array([rho_1]), os.path.join(approx_path, file))
    
    # compute relative density by monte carlo integration
    rho_2, change, n_total = monte_carlo_density(node_tags, node_coords, 
                                                 el_tags, conn, radius, x, y, 
                                                 z, dx, dy, dz, tol, 
                                                 min_points, max_points)
    print("change =", change, "n_total =", n_total)
    density_monte_carlo.append(rho_2)
    file = "montecarlo_{}.csv".format(ar)
    export_csv_file(np.array([rho_2]), os.path.join(monte_carlo_path, file))

finalize_model()

# make result plot
import matplotlib.pyplot as plt
fig, ax = plt.subplots(dpi=600)
ax.plot(aspect_ratios, density_approx, linestyle="-", linewidth=2, 
        color="#324379", label="no intersection")
ax.plot(aspect_ratios, density_monte_carlo, linestyle="-", linewidth=2, 
        color="#bf3c3c", label="Monte Carlo integration")
ax.hlines(1.0, 0.0, 1.05, linestyle="--", color="grey", linewidth=1)
ax.set_xlim([0.0, 1.05])
ax.set_ylim([0.0, 1.1])
plt.xlabel("aspect ratio")
plt.ylabel("rel. density")
plt.legend()
fig.tight_layout(pad=0.2)
plt.savefig(os.path.join(workdir, "density.pdf"))