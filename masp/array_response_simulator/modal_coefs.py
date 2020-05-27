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
#   @file   modal_coefs.py
#   @author Andrés Pérez-López
#   @date   22/08/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
from scipy.special import jv, jvp, hankel2, h2vp

from masp.array_response_simulator.sph_functions import *
from masp.validate_data_types import _validate_int, _validate_ndarray_1D, _validate_string, _validate_float


def sph_modal_coefs(N, kr, arrayType, dirCoef=None):
    """
    Modal coefficients for rigid or open spherical array

    Parameters
    ----------
    N : int
        Maximum spherical harmonic expansion order.
    kr: ndarray
        Wavenumber-radius product. Dimension = (l).
    arrayType: str
        'open', 'rigid' or 'directional'.
    dirCoef: float, optional
        Directivity coefficient of the sensor. Default to None.

    Returns
    -------
    b_N : ndarray
        Modal coefficients. Dimension = (l, N+1)

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

    """

    _validate_int('N', N, positive=True)
    _validate_ndarray_1D('kr', kr, positive=True)
    _validate_string('arrayType', arrayType, choices=['open', 'rigid', 'directional'])
    if arrayType is 'directional':
        if dirCoef is None:
            raise ValueError('dirCoef must be defined in the directional case.')
        _validate_float('dirCoef', dirCoef)

    b_N = np.zeros((kr.size, N+1), dtype='complex')

    for n in range(N+1):

        if arrayType is 'open':
            b_N[:, n] = 4 * np.pi * np.power(1j,n) * spherical_jn(n, kr)

        elif arrayType is 'rigid':
            jn = spherical_jn(n, kr)
            jnprime = spherical_jn(n, kr, derivative=True)
            hn = sph_hankel2(n, kr)
            hnprime = dsph_hankel2(n, kr)

            temp = 4 * np.pi * np.power(1j, n) * (jn - (jnprime / hnprime) * hn)
            temp[np.where(kr==0)] = 4*np.pi if n==0 else 0.
            b_N[:, n] = temp

        elif arrayType is 'directional':
            jn = spherical_jn(n, kr)
            jnprime = spherical_jn(n, kr, derivative=True)

            temp = 4 * np.pi * np.power(1j, n) * (dirCoef * jn - 1j * (1-dirCoef) * jnprime)
            b_N[:, n] = temp

    return b_N


def cyl_modal_coefs(N, kr, arrayType):
    """
    Modal coefficients for rigid or open cylindrical array

    Parameters
    ----------
    N : int
        Maximum spherical harmonic expansion order.
    kr: ndarray
        Wavenumber-radius product. Dimension = (l).
    arrayType: str
        'open' or 'rigid'.

    Returns
    -------
    b_N : ndarray
        Modal coefficients. Dimension = (l, N+1)

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    The `arrayType` options are:
    - 'open' for open array of omnidirectional sensors
    - 'rigid' for sensors mounted on a rigid baffle.

    """

    _validate_int('N', N)
    _validate_ndarray_1D('kr', kr, positive=True)
    _validate_string('arrayType', arrayType, choices=['open', 'rigid'])

    b_N = np.zeros((kr.size, N+1), dtype='complex')

    for n in range(N+1):

        if arrayType is 'open':
            b_N[:, n] = np.power(1j,n) * jv(n, kr)

        elif arrayType is 'rigid':
            jn = jv(n, kr)
            jnprime = jvp(n, kr, 1)
            hn = hankel2(n, kr)
            hnprime = h2vp(n, kr, 1)

            temp = np.power(1j, n) * (jn - (jnprime / hnprime) * hn)
            temp[np.where(kr==0)] = 1 if n==0 else 0.
            b_N[:, n] = temp

    return b_N
