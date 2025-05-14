"""
Collection of utility functions for the abqhom package.

"""

import os
import numpy as np


def voigt_to_tensor(voigt_tensor, quantity):
    """
    Transform a tensor from Voigt notation to symbolic tensor notation.

    Parameters
    ----------
    voigt_tensor : numpy.ndarray[float]
        The tensor in Voigt notation.
    quantity : str
        The quantity that is described by the tensor. Currently supported are
        'STRESS' and 'STRAIN'.

    Returns
    -------
    tensor : numpy.ndarray[float]
        The tensor in symbolic tensor notation.

    """
    quantity = quantity.upper()
    if quantity not in {"STRESS", "STRAIN"}:
        raise ValueError("Supported quantities are 'STRESS' and 'STRAIN'.")

    if voigt_tensor.shape == (6,):
        xx, yy, zz = voigt_tensor[0:3]
        yz, xz, xy = voigt_tensor[3::]

        if quantity == "STRAIN":
            yz, xz, xy = 0.5*yz, 0.5*xz, 0.5*xy

        tensor = np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])
        return tensor

    elif voigt_tensor.shape == (3,):
        xx, yy = voigt_tensor[0:2]
        xy = tensor[2]

        if quantity == "STRAIN":
            xy = 0.5*xy

        tensor = np.array([[xx, xy], [xy, yy]])
        return tensor

    else:
        raise ValueError("The shape of the provided tensor is invalid.")


def tensor_to_voigt(tensor, quantity):
    """
    Transform a tensor from symbolic tensor notation to Voigt notation.

    Parameters
    ----------
    tensor : numpy.ndarray[float]
        The tensor in symbolic tensor notation.
    quantity : str
        The quantity that is described by the tensor. Currently supported are
        'STRESS' and 'STRAIN'.

    Returns
    -------
    tensor : numpy.ndarray[float]
        The tensor in Voigt notation.

    """
    quantity = quantity.upper()
    if quantity not in {"STRESS", "STRAIN"}:
        raise ValueError("Supported quantities are 'STRESS' and 'STRAIN'.")

    if not np.allclose(tensor, tensor.T):
        raise ValueError("The provided tensor is not symmetric.")

    if tensor.shape == (3, 3):
        xx, yy, zz = tensor[0, 0], tensor[1, 1], tensor[2, 2]
        yz, xz, xy = tensor[1, 2], tensor[0, 2], tensor[0, 1]

        voigt_tensor = np.array([xx, yy, zz, yz, xz, xy])
        if quantity == "STRAIN":
            voigt_tensor[3::] *= 2
        return voigt_tensor

    elif tensor.shape == (2, 2):
        xx, yy = tensor[0, 0], tensor[1, 1]
        xy = tensor[0, 1]

        voigt_tensor = np.array([xx, yy, xy])
        if quantity == "STRAIN":
            voigt_tensor[2] *= 2
        return voigt_tensor

    else:
        raise ValueError("The shape of the provided tensor is invalid.")


def reuss(volume_fraction, C1, C2):
    """
    Compute the homogenized stiffness tensor following the Reuss theory.

    Parameters
    ----------
    volume_fraction : float
        Volume fraction of material 1.
    C1 : numpy.ndarray[float]
        Constitutive tensor of material 1 in Voigt notation.
    C2 : numpy.ndarray[float]
        Constitutive tensor of material 2 in Voigt notation.

    Returns
    -------
    C : numpy.ndarray[float]
        The homogenized constitutive tensor.

    """
    C_inv = (volume_fraction*np.linalg.inv(C1) +
             (1 - volume_fraction)*np.linalg.inv(C2))
    C = np.linalg.inv(C_inv)
    return C


def voigt(volume_fraction, C1, C2):
    """
    Compute the homogenized stiffness tensor following the Voigt theory.

    Parameters
    ----------
    volume_fraction : float
        Volume fraction of material 1.
    C1 : numpy.ndarray[float]
        Constitutive tensor of material 1 in Voigt notation.
    C2 : numpy.ndarray[float]
        Constitutive tensor of material 2 in Voigt notation.

    Returns
    -------
    C : numpy.ndarray[float]
        The homogenized constitutive tensor.

    """
    C = volume_fraction*C1 + (1 - volume_fraction)*C2
    return C


def export_to_csv_file(C, path):
    """
    Dump an array (usually a stiffness tensor) into a csv file.

    Parameters
    ----------
    C : numpy.ndarray[float]
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
    array : numpy.ndarray[dtype]
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
    files : list[str]
        List of the files found in the provided location.

    """
    if extension is None:
        files = [os.path.join(path, file) for file in os.listdir(path)]
        return files

    files = [os.path.join(path, file) for file in os.listdir(path)
             if file.endswith(extension)]
    return files
