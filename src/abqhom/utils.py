import os
import numpy as np

def reuss(vol1, C1, C2):
    """
    Compute the homogenized stiffness tensor after Reuss.

    Parameters
    ----------
    vol1 : float
        Volume fraction of material 1.
    C1 : numpy.ndarray
        Constitutive tensor of material 1 in Voigt notation.
    C2 : numpy.ndarray
        Constitutive tensor of material 2 in Voigt notation.

    Returns
    -------
    C : numpy.ndarray
        The homogenized constitutive tensor.

    """
    C_inv = vol1*np.linalg.inv(C1) + (1 - vol1)*np.linalg.inv(C2)
    C = np.linalg.inv(C_inv)
    return C

def voigt(vol1, C1, C2):
    """
    Compute the homogenized stiffness tensor after Voigt.

    Parameters
    ----------
    vol1 : float
        Volume fraction of material 1.
    C1 : numpy.ndarray
        Constitutive tensor of material 1 in Voigt notation.
    C2 : numpy.ndarray
        Constitutive tensor of material 2 in Voigt notation.

    Returns
    -------
    C : numpy.ndarray
        The homogenized constitutive tensor.

    """
    C = vol1*C1 + (1 - vol1)*C2
    return C

def export_material_tensor(C, path):
    pass

def average_stress_2d(odb, edge_length, boundary_sets, instance_name,
                      step_name, frame=-1):
    # initialize averaged stress tensor
    sig = np.zeros((2,2), dtype=float)
    
    # get internal forces
    f1_output = odb.steps[step_name].frames[frame].fieldOutputs["NFORC1"]
    f2_output = odb.steps[step_name].frames[frame].fieldOutputs["NFORC2"]
    
    # loop over all boundary nodes
    for set_name in boundary_sets:
        region = odb.rootAssembly.instances[instance_name].nodeSets[set_name]
        f1_sub = f1_output.getSubset(region=region).values
        f2_sub = f2_output.getSubset(region=region).values
        
        # collect force values
        f1_dict = {}
        f2_dict = {}
        for f1, f2 in zip(f1_sub, f2_sub):
            f1_dict[f1.nodeLabel] = f1_dict.get(f1.nodeLabel, 0.0) - f1.data
            f2_dict[f2.nodeLabel] = f2_dict.get(f2.nodeLabel, 0.0) - f2.data
            
        # update stress tensor
        for node in region.nodes:
            x = node.coordinates[0:2]
            f_inner = np.array([f1_dict[node.label], f2_dict[node.label]])
            sig += np.outer(f_inner, x)
            
    # divide by RVE volume
    sig = np.array([sig[0,0], sig[1,1], sig[0,1]])/edge_length**2
    
    return sig

# def RVE_to_matlab_input(node_tags, node_coords, el_tags, conn, file):
#     # initialize file
#     file, extension = os.path.splitext(file)
#     file += ".csv"
#     path = os.path.abspath(file)
    
#     with open(path, "w") as content:
#         # add nodes
#         content.write("//Grid,ID,x,y,z\n")
#         for tag, coords in zip(node_tags, node_coords):
#             x, y, z = coords
#             content.write("GRID,{},{},{},{}\n".format(tag, x, y, z))
        
#         # add struts
#         content.write("//Strut,ID,Start,End\n")
#         for tag, nodes in zip(el_tags, conn):
#             content.write("STRUT,{},{},{}\n".format(tag, nodes[0], nodes[1]))
            
#         return path

def read_csv_file(path, dtype=None):
    array = np.loadtxt(path, delimiter=",", dtype=dtype)
    return array

def find_all_files(path, extension=None):
    if extension is None:
        files = [os.path.join(path, file) for file in os.listdir(path)]
        return files
    
    files = [os.path.join(path, file) for file in os.listdir(path) 
             if file.endswith(extension)]
    return files