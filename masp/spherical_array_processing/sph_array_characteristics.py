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
from masp.utils import c
from masp.array_response_simulator import sph_modal_coefs
from masp.validate_data_types import _validate_float, _validate_int, _validate_string, _validate_ndarray_1D


def sph_array_noise(R, nMic, maxN, arrayType, f):
    """
    Returns noise amplification curves of a spherical mic. array

    Parameters
    ----------
    R : float
        Microphone array radius, in meter.
    nMic: int
        Number of microphone capsules.
    maxN : int
        Maximum spherical harmonic expansion order.
    arrayType: str
        'open' or 'rigid'.
    f: ndarray
        Frequencies where to perform estimation. Dimension = (l).

    Returns
    -------
    g2 : ndarray
        Noise amplification curve. Dimension = (l, maxN)
    g2_lin: ndarray
        Linear log-log interpolation of noise. Dimension = (l, maxN)

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
