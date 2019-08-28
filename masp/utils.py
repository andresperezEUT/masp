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
#   @file   utils.py
#   @author Andrés Pérez-López
#   @date   30/07/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import numpy as np
from scipy.special import sph_harm
from math import factorial

from masp.validate_data_types import _validate_int, _validate_ndarray_2D, _validate_string, _validate_ndarray_1D

c = 343.
C = 3

def cart2sph(x, y, z):
    """
    TODO
    implemented from matlab
    :param x:
    :param y:
    :param z:
    :return:
    """
    hypotxy = np.hypot(x, y)
    r = np.hypot(hypotxy, z)
    elev = np.arctan2(z, hypotxy)
    az = np.arctan2(y, x)
    return az, elev, r

def sph2cart(az,elev,r):
    """
    TODO
    implemented from matlab
    :param az:
    :param elev:
    :param r:
    :return:
    """
    z = r * np.sin(elev)
    rcoselev = r * np.cos(elev)
    x = rcoselev * np.cos(az)
    y = rcoselev * np.sin(az)
    return np.asarray([x, y, z]).T


def get_sh(N, dirs, basisType):
    """
    Get spherical harmonics up to order N evaluated at given angular positions.

    Parameters
    ----------
    N : int
        Maximum spherical harmonic expansion order.
    dirs : ndarray
        Evaluation directions. Dimension = (nDirs, 2).
        Directions are expected in radians, expressed in pairs [azimuth, inclination].

    basisType : str
        Type of spherical harmonics. Either 'complex' or 'real'.

    Returns
    -------
    Y : ndarray
        Wall absorption coefficients . Dimension = (nDirs, nHarm).

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.
    NotImplementedError: for basisType = 'complex'

    Notes
    -----
    Ouput dimension is given by: nHarm = (N+1)^2

    Inclination is defined as the angle from zenith: inclination = pi/2-elevation

    TODO: implement complex basis?
    """

    _validate_int('order', N, positive=True)
    _validate_ndarray_2D('dirs', dirs, shape1=2)
    _validate_string('basisType', basisType, choices=['complex', 'real'])

    nDirs = dirs.shape[0]
    nHarm = np.power(N+1, 2)
    Y_N = np.zeros((nDirs, nHarm))

    def delta_kronecker(q1, q2):
        return 1 if q1 == q2 else 0

    def norm(m):
        """
        TODO
        SN3D FACTOR, REMOVE CONDON SHOTLEY
        IT MUST BE MULTIPLIED BY sqrt(4PI) to be normalized to 1
        :param m:
        :return:
        """
        return np.power(-1, np.abs(m)) * np.sqrt(2 - delta_kronecker(0, np.abs(m)))

    if basisType is 'complex':
        raise NotImplementedError

    elif basisType is 'real':
        harm_idx = 0
        for l in range(N + 1):
            for m in range(-l,0):
                Y_N[:, harm_idx] = np.imag(sph_harm(np.abs(m), l, dirs[:, 0], dirs[:, 1])) * norm(m)
                harm_idx += 1
            for m in range(l + 1):
                Y_N[:, harm_idx] = np.real(sph_harm(m, l, dirs[:, 0], dirs[:, 1])) * norm(m)
                harm_idx += 1

    return Y_N


def lagrange(N, delays):
    """
    Design a fractional delay order-N filter matrix with polynomial interpolation.

    Parameters
    ----------
    N : int
        Filter order.
    delays : ndarray
       Target fractional delays, in samples. Dimension = (1).

    Returns
    -------
    h : ndarray
        Target filter. Dimension = (N+1, len(delays))

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    For best results, delay should be near N/2 +/- 1.
    """

    _validate_int('N', N, positive=True)
    _validate_ndarray_1D('delays', delays, positive=True)

    n = np.arange(N+1)
    h = np.ones((N+1, delays.size))
    for l in range(delays.size):
        for k in range(N+1):
            idx = n[n != k]
            h[idx, l] = h[idx, l] * (delays[l]-k) / (n[idx]-k)
    return h
