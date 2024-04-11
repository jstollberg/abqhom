import os
import sys

# import gmsh and homogenization package
base_path = "C:/Users/jonat/Documents/"
gmsh_path = os.path.join(base_path, "gmsh", "lib")
abqhom_path = os.path.join(base_path, 
                            "Institute for Mechanics", 
                            "abqhom", 
                            "src")
sys.path.append(gmsh_path)
sys.path.append(abqhom_path)
from abqhom.RVE import write_matlab_input, finalize_model
from abqhom.examples import simple_RVE_3d
import gmsh

dx = dy = dz = 10
lc = 5

model_name, group_map = simple_RVE_3d(dx=dx, dy=dy, dz=dz, lc=lc)
gmsh.fltk.run()

write_matlab_input(model_name, "stochastic.txt")
finalize_model()