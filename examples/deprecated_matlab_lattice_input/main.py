import os
import sys

# generate input file for the matlab homogenization tool available from:
# https://www.mathworks.com/matlabcentral/fileexchange/67457-3d-homogenization-of-cellular-materials

base_path = "C:/Users/jonat/Documents/"
gmsh_path = os.path.join(base_path, "gmsh", "lib")
abqhom_path = os.path.join(base_path,
                           "abqhom",
                           "src")
sys.path.append(gmsh_path)
sys.path.append(abqhom_path)
from abqhom.RVE import write_matlab_input, finalize_model
from abqhom.examples import simple_RVE_3d
import gmsh

# RVE size and strut length
dx = dy = dz = 10
lc = 5

# write input file
model_name, group_map = simple_RVE_3d(dx=dx, dy=dy, dz=dz, lc=lc)
write_matlab_input(model_name, "stochastic.txt")
# gmsh.fltk.run()
finalize_model()