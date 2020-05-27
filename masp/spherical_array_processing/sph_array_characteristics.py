# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Copyright (c) 2019, Eurecat / UPF
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the <organization> nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
#   @file   sph_array_characteristics.py
#   @author Andrés Pérez-López
#   @date   27/09/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import numpy as np
from masp.utils import c, C, check_cond_number_sht, elev2incl
from masp.array_response_simulator import sph_modal_coefs
from masp.validate_data_types import _validate_float, _validate_int, _validate_string, _validate_ndarray_1D, \
    _validate_ndarray_2D


def sph_array_noise(R, nMic, maxN, arrayType, f):
    """
    Returns noise amplification curves of a spherical mic. array

    Parameters
    ----------
    R : float
        Microphone array radius, in meter.
    nMic : int
        Number of microphone capsules.
    maxN : int
        Maximum spherical harmonic expansion order.
    arrayType : str
        'open' or 'rigid'.
    f : ndarray
        Frequencies where to perform estimation. Dimension = (l).

    Returns
    -------
    g2 : ndarray
        Noise amplification curve. Dimension = (l, maxN)
    g2_lin: ndarray
        Linear log-log interpolation of noise. Dimension = (l, maxN).

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    The `arrayType` options are:
    - 'open' for open array of omnidirectional sensors, or
    - 'rigid' for sensors mounted on a rigid baffle.

    `g2_lin` is an approximation of the curves at low frequencies showing
    the linear behaviour in the log-log axis, with a 6n dB slope.

    """
    
    _validate_float('R', R, positive=True)
    _validate_int('nMic', nMic, positive=True)
    _validate_int('maxN', maxN, positive=True)
    _validate_string('arrayType', arrayType, choices=['open', 'rigid'])
    _validate_ndarray_1D('f', f, positive=True)

    # Frequency axis
    kR = 2 * np.pi * f * R / c
    # Modal responses
    bN = sph_modal_coefs(maxN, kR, arrayType) / (4*np.pi)
    # Noise power response
    g2 = 1. / ( nMic * np.power(np.abs(bN), 2) )

    # Approximate linearly
    p = -(6 / 10) * np.arange(1, maxN+1) / np.log10(2)
    bN_lim0 = (sph_modal_coefs(maxN, np.asarray([1]), arrayType) / (4 * np.pi)).squeeze()
    a = 1. / ( nMic * np.power(np.abs(bN_lim0[1:]), 2) )

    g2_lin = np.zeros((kR.size, maxN))
    for n in range(maxN):
        g2_lin[:, n] = a[n] * np.power(kR, p[n])

    return g2, g2_lin


def sph_array_noise_threshold(R, nMic, maxG_db, maxN, arrayType, dirCoef=None):
    """
    Returns frequency limits for noise amplification of a spherical mic. array

    Parameters
    ----------
    R : float
        Microphone array radius, in meter.
    nMic : int
        Number of microphone capsules.
    maxG_db : float
        max allowed amplification for the noise level. maxG_db = 20*log10(maxG)
    maxN : int
        Maximum spherical harmonic expansion order.
    arrayType: str
        'open', 'rigid' or 'directional'
    dirCoef: float, optional
        Directivity coefficient of the sensor. Default to None.

    Returns
    -------
    f_lim : ndarray
       Frequency points where threhsold is reached, for orders 1:maxN. Dimension = (maxN).

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    The `arrayType` options are:
    - 'open' for open array of omnidirectional sensors,
    - 'rigid' for sensors mounted on a rigid baffle,
    - 'directional' for an open array of first-order directional microphones determined by `dirCoef`.

    `dirCoef` is relevant (and required) only in the 'directional' type.
    `dirCoef` ranges from 0 (omni) to 1 (dipole), where for example 0.5 is a cardioid sensor.
    In the 0 case it is equivalent to an open array of omnis.
    The first order directivity function is defined as d(theta) = dirCoeff + (1-dirCoeff)*cos(theta).


    The method returns the frequencies that the noise in the
    output channels of a SMA, after performing the SHT and equalization of
    the output signals, reaches a certain user-defined threshold maxG_db.
    The frequencies are computed only at the lower range of each order,
    where its response decays rapidly, ignoring for example the nulls of an
    open array at the higher frequencies. The estimation of the limits are
    based on a linear approximation of the log-log response found e.g. in

        Sector-based Parametric Sound Field Reproduction in the Spherical Harmonic Domain
        A Politis, J Vilkamo, V Pulkki
        IEEE Journal of Selected Topics in Signal Processing 9 (5), 852 - 866

    """

    _validate_float('R', R, positive=True)
    _validate_int('nMic', nMic, positive=True)
    _validate_float('maxG_db', maxG_db)
    _validate_int('maxN', maxN, positive=True)
    _validate_string('arrayType', arrayType, choices=['open', 'rigid', 'directional'])
    if arrayType is 'directional':
        if dirCoef is None:
            raise ValueError('dirCoef must be defined in the directional case.')
        _validate_float('dirCoef', dirCoef)

    f_lim = np.zeros(maxN)
    for n in range(1, maxN+1):
        bn = sph_modal_coefs(n, np.asarray([1]), arrayType, dirCoef) / (4 * np.pi)
        bn = bn.flatten()[-1]
        maxG = np.power(10, (maxG_db / 10))
        kR_lim = np.power( (maxG * nMic * np.power(np.abs(bn), 2)), (-10 * np.log10(2) / (6 * n)) )
        f_lim[n-1] = kR_lim * c / (2 * np.pi * R)

    return f_lim


def sph_array_alias_lim(R, nMic, maxN, mic_dirs_rad, mic_weights=None):
    """
    Get estimates of the aliasing limit of a spherical array, in three different ways.

    Parameters
    ----------
    R : float
        Microphone array radius, in meter.
    nMic : int
        Number of microphone capsules.
    maxN : int
        Maximum spherical harmonic expansion order.
    mic_dirs_rad : ndarray
        Evaluation directions. Dimension = (nDirs, 2).
        Directions are expected in radians, expressed in pairs [azimuth, elevation].
    dirCoef: ndarray, optional
       Vector of weights used to improve orthogonality of the SH transform. Dimension = (nDirs)

    Returns
    -------
    f_alias : ndarray
       the alising limit estimates. Dimension = (3).

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    First estimate takes into account only the radius and a nominal order
    that the array is expected to support, it is the simplest one and it
    expresses the kR = maxN rule.
    The second estimate is based on the number of microphones, and it can
    be more relaxed than the first, if the nominal supported order is less
    than maxN<floor(sqrt(Nmic)-1).
    The third takes into account microphone numbers and directions, and it
    is based on the orthogonality of the SH matrix for the microphone
    positions expressed through the condition number.

    # todo check RETURN VALUES (need for very big rtol if cond_N returned)
    """

    _validate_float('R', R, positive=True)
    _validate_int('nMic', nMic, positive=True)
    _validate_int('maxN', maxN, positive=True)
    _validate_ndarray_2D('mic_dirs_rad', mic_dirs_rad, shape1=C-1)
    if mic_weights is not None:
        _validate_ndarray_1D('mic_weights', mic_weights, size=mic_dirs_rad.shape[0])

    f_alias = np.zeros(3)

    # Conventional kR = N assumption
    f_alias[0] = c * maxN / (2 * np.pi * R)

    # Based on the floor of the number of microphones, uniform arrangement
    f_alias[1] = c * np.floor(np.sqrt(nMic) - 1) / (2 * np.pi * R)

    # Based on condition number of the SHT matrix
    maxOrder = int(np.ceil(np.sqrt(nMic) - 1))

    cond_N = check_cond_number_sht(maxOrder, elev2incl(mic_dirs_rad), 'real', mic_weights)
    trueMaxOrder = np.flatnonzero(cond_N<np.power(10,4))[-1]  # biggest element passing the condition
    f_alias[2] = c * trueMaxOrder / (2 * np.pi * R)

    return f_alias


