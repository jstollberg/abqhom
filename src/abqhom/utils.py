import os
import numpy as np

def voigt_to_strain_tensor(strain):
    """
    Convert strain tensor from Voigt notation to tensor notation.

    Parameters
    ----------
    strain : numpy.ndarray
        Strain tensor in Voigt notation.

    Returns
    -------
    eps : numpy.ndarray
        Strain tensor in symbolic notation.

    """
    if strain.shape == (6,):
        eps_xx, eps_yy, eps_zz = strain[0:3]
        eps_yz, eps_xz, eps_xy = strain[3::]/2
        
        eps = np.array([[eps_xx, eps_xy, eps_xz],
                        [eps_xy, eps_yy, eps_yz],
                        [eps_xz, eps_yz, eps_zz]])
        
    elif strain.shape == (3,):
        eps_xx, eps_yy = strain[0:2]
        eps_xy = strain[2]/2
        
        eps = np.array([[eps_xx, eps_xy],
                        [eps_xy, eps_yy]])
        
    else:
        raise ValueError("strain must be provided in Voigt notation")
        
    return eps

def voigt_to_stress_tensor(stress):
    """
    Convert stress tensor in Voigt notation to tensor notation.

    Parameters
    ----------
    stress : numpy.ndarray
        Stress tensor in Voigt notation.

    Returns
    -------
    sig : numpy.ndarray
        Stress tensor in symbolic notation.

    """
    if stress.shape == (6,):
        sig_xx, sig_yy, sig_zz, sig_yz, sig_xz, sig_xy = stress
        
        sig = np.array([[sig_xx, sig_xy, sig_xz],
                        [sig_xy, sig_yy, sig_yz],
                        [sig_xz, sig_yz, sig_zz]])
        
    elif stress.shape == (3,):
        sig_xx, sig_yy, sig_xy = stress
        
        sig = np.array([[sig_xx, sig_xy],
                        [sig_xy, sig_yy]])
        
    else:
        raise ValueError("stress must be provided in Voigt notation")
        
    return sig
        
def strain_tensor_to_voigt(strain):
    """
    Convert strain tensor from tensor notation to Voigt notation.

    Parameters
    ----------
    strain : numpy.ndarray
        Strain tensor in symbolic notation.

    Returns
    -------
    eps : numpy.ndarray
        Strain tensor in Voigt notation.

    """
    if strain.shape == (3,3):
        eps_xx, eps_yy, eps_zz = strain[0,0], strain[1,1], strain[2,2]
        eps_yz, eps_xz, eps_xy = strain[1,2], strain[0,2], strain[0,1]
        
        eps = np.array([eps_xx, eps_yy, eps_zz, eps_yz, eps_xz, eps_xy])
        eps[3::] *= 0.5
        
    elif strain.shape == (2,2):
        eps_xx, eps_yy = strain[0,0], strain[1,1]
        eps_xy = strain[0,1]
        
        eps = np.array([eps_xx, eps_yy, 0.5*eps_xy])
        
    else:
        raise ValueError("input has wrong shape")
        
    return eps
        
def stress_tensor_to_voigt(stress):
    """
    Convert stress tensor from tensor notation to Voigt notation.

    Parameters
    ----------
    stress : numpy.ndarray
        Stress tensor in symbolic notation.

    Returns
    -------
    sig : numpy.ndarray
        Stress tensor in Voigt notation.

    """
    if stress.shape == (3,3):
        sig_xx, sig_yy, sig_zz = stress[0,0], stress[1,1], stress[2,2]
        sig_yz, sig_xz, sig_xy = stress[1,2], stress[0,2], stress[0,1]
        
        sig = np.array([sig_xx, sig_yy, sig_zz, sig_yz, sig_xz, sig_xy])
        
    elif stress.shape == (2,2):
        sig_xx, sig_yy = stress[0,0], stress[1,1]
        sig_xy = stress[0,1]
        
        sig = np.array([sig_xx, sig_yy, sig_xy])
        
    else:
        raise ValueError("input has wrong shape")
        
    return sig

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

def export_csv_file(C, path):
    """
    Dump an array into a csv file.

    Parameters
    ----------
    C : numpy.ndarray
        The array to export.
    path : str
        The file path.

    Returns
    -------
    None.

    """
    path, extension = os.path.splitext(path)
    path += ".csv"
    np.savetxt(path, C, delimiter=",")
    
def read_csv_file(path, dtype=None):
    """
    Read a csv file into an array.

    Parameters
    ----------
    path : str
        Path to the file.
    dtype : type, optional
        Data type of the array. The default is None.

    Returns
    -------
    array : numpy.ndarray
        The array containing the values from the csv file.

    """
    array = np.loadtxt(path, delimiter=",", dtype=dtype)
    return array

def find_all_files(path, extension=None):
    """
    Find all files in a specific location.

    Parameters
    ----------
    path : str
        The location to browse.
    extension : str, optional
        The file type to look for. The default is None.

    Returns
    -------
    files : list
        List of the files found in the provided location.

    """
    if extension is None:
        files = [os.path.join(path, file) for file in os.listdir(path)]
        return files
    
    files = [os.path.join(path, file) for file in os.listdir(path) 
             if file.endswith(extension)]
    return files

def best_isotropic_approximation(C):
    # if C.shape == (3,3):
    #     C = np.array([[C[0,0], C[0,1], ]])
    
    lam0 = 1/15*(C[0,0] + C[1,1] + C[2,2] - 2*(C[3,3] + C[4,4] + C[5,5])
                 + 4*(C[0,1] + C[0,2] + C[1,2]))
    mu0 = 1/15*(C[0,0] + C[1,1] + C[2,2] + 3*(C[3,3] + C[4,4] + C[5,5])
                - (C[0,1] + C[0,2] + C[1,2]))

    # C_iso = lambda lam, mu: np.array([[lam + 2*mu, lam, lam, 0.0, 0.0, 0.0],
    #                                   [lam, lam + 2*mu, lam, 0.0, 0.0, 0.0],
    #                                   [lam, lam, lam + 2*mu, 0.0, 0.0, 0.0],
    #                                   [0.0, 0.0, 0.0, mu, 0.0, 0.0],
    #                                   [0.0, 0.0, 0.0, 0.0, mu, 0.0],
    #                                   [0.0, 0.0, 0.0, 0.0, 0.0, mu]])

def anisotropy_index():
    pass



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