# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Collection of utility functions for the abqhom package.

"""

import os
import numpy as np


def voigt_to_tensor(voigt_tensor, quantity):
    r"""
    Convert symmetric stress or strain data from Voigt to tensor form.

    The three-dimensional ordering is
    :math:`(11,22,33,23,13,12)` and the two-dimensional ordering is
    :math:`(11,22,12)`. Strain input uses engineering shear components,
    while stress input uses tensorial shear components.

    Parameters
    ----------
    voigt_tensor : ndarray
        One Voigt vector or an array of Voigt vectors.
    quantity : str
        The quantity that is described by the tensor. Currently
        supported are 'STRESS' and 'STRAIN'.

    Returns
    -------
    tensor : ndarray
        One symmetric second-order tensor or an array of such tensors.

    """
    quantity = quantity.upper()
    if quantity not in {"STRESS", "STRAIN"}:
        raise ValueError("Supported quantities are 'STRESS' and 'STRAIN'.")

    if len(voigt_tensor.shape) == 1:
        voigt_tensor = np.expand_dims(voigt_tensor, 0)

    tensor = []
    for t in voigt_tensor:
        if t.shape == (6,):
            xx, yy, zz = t[0:3]
            yz, xz, xy = t[3::]

            if quantity == "STRAIN":
                yz, xz, xy = 0.5*yz, 0.5*xz, 0.5*xy

            tensor.append(np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]]))

        elif t.shape == (3,):
            xx, yy = t[0:2]
            xy = t[2]

            if quantity == "STRAIN":
                xy = 0.5*xy

            tensor.append(np.array([[xx, xy], [xy, yy]]))

        else:
            raise ValueError("The shape of the provided tensor is invalid.")

    tensor = np.array(tensor)
    if tensor.shape[0] == 1:
        tensor = tensor[0]
    return tensor


def tensor_to_voigt(tensor, quantity, check_symmetry=True):
    r"""
    Convert symmetric stress or strain data from tensor to Voigt form.

    The three-dimensional ordering is
    :math:`(11,22,33,23,13,12)` and the two-dimensional ordering is
    :math:`(11,22,12)`. Engineering shear components are returned for
    strain, while stress retains tensorial shear components.

    Parameters
    ----------
    tensor : ndarray
        One second-order tensor or an array of such tensors.
    quantity : str
        The quantity that is described by the tensor. Currently
        supported are 'STRESS' and 'STRAIN'.
    check_symmetry : bool
        If True, the input will be checked for symmetry. The default
        is True.

    Returns
    -------
    voigt_tensor : ndarray
        One Voigt vector or an array of Voigt vectors.

    """
    quantity = quantity.upper()
    if quantity not in {"STRESS", "STRAIN"}:
        raise ValueError("Supported quantities are 'STRESS' and 'STRAIN'.")

    if len(tensor.shape) == 2:
        tensor = np.expand_dims(tensor, 0)

    voigt = []
    for t in tensor:
        if not np.allclose(t, t.T) and check_symmetry:
            raise ValueError("The provided tensor is not symmetric.")

        if t.shape == (3, 3):
            xx, yy, zz = t[0, 0], t[1, 1], t[2, 2]
            yz, xz, xy = t[1, 2], t[0, 2], t[0, 1]

            voigt_tensor = np.array([xx, yy, zz, yz, xz, xy])
            if quantity == "STRAIN":
                voigt_tensor[3::] *= 2
            voigt.append(voigt_tensor)

        elif t.shape == (2, 2):
            xx, yy = t[0, 0], t[1, 1]
            xy = t[0, 1]

            voigt_tensor = np.array([xx, yy, xy])
            if quantity == "STRAIN":
                voigt_tensor[2] *= 2
            voigt.append(voigt_tensor)

        else:
            raise ValueError("The shape of the provided tensor is invalid.")

    voigt = np.array(voigt, dtype=float)
    if voigt.shape[0] == 1:
        voigt = voigt[0]

    return voigt


def hill_mandel_residual(macroscopic_power, microscopic_power,
                         power_tolerance):
    r"""
    Compute the normalized discrete Hill--Mandel residual.

    For the macroscopic and microscopic incremental-work densities
    :math:`\Delta w_\mathrm{macro}` and
    :math:`\Delta w_\mathrm{micro}`, respectively, the residual is

    .. math::
        \epsilon_\mathrm{HM} =
        \frac{|\Delta w_\mathrm{macro}-\Delta w_\mathrm{micro}|}
        {\max(|\Delta w_\mathrm{macro}|,
              |\Delta w_\mathrm{micro}|,\Delta w_\mathrm{tol})}.

    Parameters
    ----------
    macroscopic_power : float or ArrayLike
        Macroscopic incremental-work density
        :math:`\Delta w_\mathrm{macro}`.
    microscopic_power : float or ArrayLike
        Microscopic incremental-work density
        :math:`\Delta w_\mathrm{micro}`.
    power_tolerance : float
        Positive regularizing incremental-work-density scale
        :math:`\Delta w_\mathrm{tol}`.

    Returns
    -------
    residual : float or ndarray
        The normalized residual. Array inputs are evaluated elementwise.

    """
    power_tolerance = float(power_tolerance)
    if power_tolerance <= 0.0:
        raise ValueError("The power tolerance must be positive.")

    macroscopic_power = np.asarray(macroscopic_power, dtype=float)
    microscopic_power = np.asarray(microscopic_power, dtype=float)
    scale = np.maximum(np.maximum(np.abs(macroscopic_power),
                                  np.abs(microscopic_power)),
                       power_tolerance)
    residual = np.abs(macroscopic_power - microscopic_power)/scale

    if residual.ndim == 0:
        return float(residual)
    return residual


def reuss(volume_fraction, C1, C2):
    r"""
    Return the Reuss estimate of the effective stiffness matrix.

    For the volume fraction :math:`c` of phase 1,
    :math:`\bar{\mathbb{C}}_{\mathrm{R}}^{\mathrm{V}} =
    (c(\mathbb{C}^{\mathrm{V}}_1)^{-1}
    +(1-c)(\mathbb{C}^{\mathrm{V}}_2)^{-1})^{-1}`.

    Parameters
    ----------
    volume_fraction : float
        Volume fraction of phase 1.
    C1 : ndarray
        Stiffness matrix of phase 1 in Voigt notation.
    C2 : ndarray
        Stiffness matrix of phase 2 in Voigt notation.

    Returns
    -------
    C : ndarray
        Reuss estimate of the effective stiffness matrix.

    """
    C_inv = (volume_fraction*np.linalg.inv(C1) +
             (1 - volume_fraction)*np.linalg.inv(C2))
    C = np.linalg.inv(C_inv)
    return C


def voigt(volume_fraction, C1, C2):
    r"""
    Return the Voigt estimate of the effective stiffness matrix.

    For the volume fraction :math:`c` of phase 1,
    :math:`\bar{\mathbb{C}}_{\mathrm{V}}^{\mathrm{V}} =
    c\mathbb{C}^{\mathrm{V}}_1+(1-c)\mathbb{C}^{\mathrm{V}}_2`.

    Parameters
    ----------
    volume_fraction : float
        Volume fraction of phase 1.
    C1 : ndarray
        Stiffness matrix of phase 1 in Voigt notation.
    C2 : ndarray
        Stiffness matrix of phase 2 in Voigt notation.

    Returns
    -------
    C : ndarray
        Voigt estimate of the effective stiffness matrix.

    """
    C = (volume_fraction*C1) + (1 - volume_fraction)*C2
    return C


def export_to_csv_file(C, path):
    """
    Write an array to a CSV file.

    Parameters
    ----------
    C : ArrayLike
        The array to export.
    path : str or os.PathLike
        The file path.

    Returns
    -------
    None

    """
    path, _ = os.path.splitext(path)
    path += ".csv"
    np.savetxt(path, C, delimiter=",")


def read_csv_file(path, dtype=None):
    """
    Read a CSV file into an array.

    Parameters
    ----------
    path : str or os.PathLike
        Path to the file.
    dtype : data-type, optional
        Any data-type specification accepted by ``numpy.loadtxt``. The
        default is None.

    Returns
    -------
    array : ndarray
        The array containing the values from the CSV file.

    """
    array = np.loadtxt(path, delimiter=",", dtype=dtype)
    return array
