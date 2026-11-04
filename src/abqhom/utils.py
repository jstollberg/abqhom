# SPDX-FileCopyrightText: 2023-2026 Jonathan Stollberg
# SPDX-License-Identifier: GPL-3.0-or-later

"""Numerical utilities for homogenization post-processing."""

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
    check_symmetry : bool, optional
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

    voigt_data = []
    for t in tensor:
        if check_symmetry and not np.allclose(t, t.T):
            raise ValueError("The provided tensor is not symmetric.")

        if t.shape == (3, 3):
            xx, yy, zz = t[0, 0], t[1, 1], t[2, 2]
            yz, xz, xy = t[1, 2], t[0, 2], t[0, 1]

            voigt_tensor = np.array([xx, yy, zz, yz, xz, xy])
            if quantity == "STRAIN":
                voigt_tensor[3::] *= 2
            voigt_data.append(voigt_tensor)

        elif t.shape == (2, 2):
            xx, yy = t[0, 0], t[1, 1]
            xy = t[0, 1]

            voigt_tensor = np.array([xx, yy, xy])
            if quantity == "STRAIN":
                voigt_tensor[2] *= 2
            voigt_data.append(voigt_tensor)

        else:
            raise ValueError("The shape of the provided tensor is invalid.")

    voigt_data = np.array(voigt_data, dtype=float)
    if voigt_data.shape[0] == 1:
        voigt_data = voigt_data[0]

    return voigt_data


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


def _validate_constituent_data(volume_fractions, stiffnesses):
    """Validate and return constituent fractions and stiffness matrices."""

    volume_fractions = np.asarray(volume_fractions, dtype=float)
    if volume_fractions.ndim != 1 or volume_fractions.size == 0:
        raise ValueError("Volume fractions must be a nonempty 1D sequence.")
    if (not np.all(np.isfinite(volume_fractions))
            or np.any(volume_fractions < 0.0)):
        raise ValueError("Volume fractions must be finite and nonnegative.")
    if not np.isclose(np.sum(volume_fractions), 1.0,
                      rtol=1.0e-12, atol=1.0e-12):
        raise ValueError("Volume fractions must add up to one.")

    stiffnesses = np.asarray(stiffnesses, dtype=float)
    if (stiffnesses.ndim != 3
            or stiffnesses.shape[0] != volume_fractions.size
            or stiffnesses.shape[1] != stiffnesses.shape[2]):
        raise ValueError(
            "Provide one square stiffness matrix per constituent.")
    if not np.all(np.isfinite(stiffnesses)):
        raise ValueError("Stiffness matrices must contain finite values.")

    return volume_fractions, stiffnesses


def reuss(volume_fractions, stiffnesses):
    r"""
    Return the Reuss estimate of the effective stiffness matrix.

    For :math:`N` constituents with volume fractions :math:`c_a`,
    :math:`\bar{\mathbb{C}}_{\mathrm{R}}^{\mathrm{V}} =
    \left(\sum_{a=1}^N c_a
    (\mathbb{C}^{\mathrm{V}}_a)^{-1}\right)^{-1}`.

    Parameters
    ----------
    volume_fractions : ArrayLike
        Nonnegative constituent volume fractions that add up to one.
    stiffnesses : ArrayLike
        One square stiffness matrix per constituent in Voigt notation.

    Returns
    -------
    C : ndarray
        Reuss estimate of the effective stiffness matrix.

    """
    volume_fractions, stiffnesses = _validate_constituent_data(
        volume_fractions, stiffnesses)
    compliance = np.zeros_like(stiffnesses[0])
    for volume_fraction, stiffness in zip(volume_fractions, stiffnesses):
        if volume_fraction > 0.0:
            compliance += volume_fraction*np.linalg.inv(stiffness)
    return np.linalg.inv(compliance)


def voigt(volume_fractions, stiffnesses):
    r"""
    Return the Voigt estimate of the effective stiffness matrix.

    For :math:`N` constituents with volume fractions :math:`c_a`,
    :math:`\bar{\mathbb{C}}_{\mathrm{V}}^{\mathrm{V}} =
    \sum_{a=1}^N c_a\mathbb{C}^{\mathrm{V}}_a`.

    Parameters
    ----------
    volume_fractions : ArrayLike
        Nonnegative constituent volume fractions that add up to one.
    stiffnesses : ArrayLike
        One square stiffness matrix per constituent in Voigt notation.

    Returns
    -------
    C : ndarray
        Voigt estimate of the effective stiffness matrix.

    """
    volume_fractions, stiffnesses = _validate_constituent_data(
        volume_fractions, stiffnesses)
    return np.tensordot(volume_fractions, stiffnesses, axes=(0, 0))


def export_to_csv_file(array, path, header=None, fmt="%.18e"):
    """
    Write an array to a CSV file.

    Parameters
    ----------
    array : ArrayLike
        The array to export.
    path : str or os.PathLike
        The file path. An existing suffix is replaced by ``.csv``.
    header : str, optional
        Comma-separated column names written as the first row. The default
        is None.
    fmt : str or sequence of str, optional
        Numeric format accepted by ``numpy.savetxt``. The default is
        ``"%.18e"``.

    """
    path, _ = os.path.splitext(path)
    path += ".csv"
    options = {"delimiter": ",", "fmt": fmt}
    if header is not None:
        options.update(header=header, comments="")
    np.savetxt(path, array, **options)


def read_csv_file(path, dtype=None, has_header=False):
    """
    Read a CSV file into an array.

    Parameters
    ----------
    path : str or os.PathLike
        Path to the file.
    dtype : data-type, optional
        Any data-type specification accepted by ``numpy.loadtxt``. The
        default is None.
    has_header : bool, optional
        If True, parse the first row as comma-separated column names and
        return it together with the array. The default is False.

    Returns
    -------
    array : ndarray
        The array containing the values from the CSV file.
    header : tuple of str
        Column names. Returned only if ``has_header`` is True.

    """
    skiprows = 0
    header = None
    if has_header:
        with open(path, "r") as csv_file:
            first_row = csv_file.readline()
        if not first_row:
            raise ValueError("A CSV header was requested from an empty file.")
        header = tuple(
            name.strip()
            for name in first_row.lstrip("# ").strip().split(","))
        skiprows = 1

    array = np.loadtxt(
        path, delimiter=",", dtype=dtype, skiprows=skiprows)
    return (array, header) if has_header else array
