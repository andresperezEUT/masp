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

from masp.validate_data_types import _validate_int, _validate_ndarray_2D, _validate_string

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

    Inclination is defined as the angle from zenit: inclination = pi/2-elevation
    """

    _validate_int('order', N, positive=True)
    _validate_ndarray_2D('dirs', dirs, shape1=2)
    _validate_string('basisType', basisType, choices=['complex', 'real'])

    # expand 1dim
    if np.ndim(dirs) == 1:
        dirs = dirs[np.newaxis, :]

    nDirs = dirs.shape[0]
    nHarm = np.power(N+1, 2)

    Y_N = np.zeros((nHarm, nDirs))

    if basisType is 'complex':
        raise NotImplementedError

    elif basisType is 'real':
        for dir_idx, dir in enumerate(dirs):
            azi = dir[0]
            incl = np.pi/2 - dir[1]
            Y_N[:, dir_idx] = get_ambisonics_coefs(azi, incl, N)
        return Y_N.transpose()


def delta_kronecker(q1, q2):
    if (q1==q2): return 1
    else:        return 0

def get_real_spherical_harmonic(azimuth, elevation, ambisonics_order, ambisonics_degree):
    # NOTE THAT EVERYTHING IS CHANGED RESPECT TO THE MAN ENTRY
    # here, we use phi as azimuth and theta as elevation
    # furthermore, L is ambisonics order and M is ambisonics degree
    return np.real(sph_harm(ambisonics_degree, ambisonics_order, azimuth, elevation))

def get_imag_spherical_harmonic(azimuth, elevation, ambisonics_order, ambisonics_degree):
    # NOTE THAT EVERYTHING IS CHANGED RESPECT TO THE MAN ENTRY
    # here, we use phi as azimuth and theta as elevation
    # furthermore, L is ambisonics order and M is ambisonics degree
    return np.imag(sph_harm(ambisonics_degree,ambisonics_order,azimuth,elevation))


def get_spherical_harmonic_normalization_coef(order,degree):
    # standard normalization function, without Condon-Shotley

    # we can use it to divide the sph_harm function by this number and get 1-normalization

    # this is the default coef given by sph_harm function

    # this is an implementation of the normalization function provided by sph_harm
    # see https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.special.sph_harm.html
    #
    # it follows the standard normalization convention as used IIRC in physics
    # since we want SN3D and not this one, we use this function to compute it
    # in order to substract it to the resulting coefficients, and then applying
    # the SN3D normalization
    #
    # notice that Condon-Shortley phase is not included here, since it is calculated internally
    # by means of lpmv function
    l = order
    m = abs(degree)
    return np.sqrt( ( ((2*l)+1) / (4*np.pi) ) * ( float(factorial(l-m)) / float(factorial(l+m)) ) )


def get_spherical_harmonic(azimuth, elevation, ambisonics_order, ambisonics_degree):

    l = ambisonics_order
    m = abs(ambisonics_degree)

    # in the sph_harm, the elevation is defined in the range [0..pi]
    # luckily the sph_harm reference system follows the right-hand rule
    # so there is no other changes involved
    elevation_internal = elevation - (np.pi/2)

    # the harmonics with positive degree are the real part of the solution
    # and the ones with negative degree are the imaginary part
    # degree 0 is symmetrical for both solutions
    # check http://mathworld.wolfram.com/SphericalHarmonic.html
    if (ambisonics_degree >= 0):
        coef = get_real_spherical_harmonic(azimuth,elevation_internal,l, m)
    else:
        coef = get_imag_spherical_harmonic(azimuth,elevation_internal,l, m)

    # normalization
    #
    # get_real_spherical_harmonic returns a value with a different normalization
    # that the one we want
    # (check https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.special.sph_harm.html)
    #
    # therefore we "take it out" from the result, and apply our desired SN3D normalization
    # as described in Daniel2003, 2.1 (http://gyronymo.free.fr/audio3D/publications/AES23%20NFC%20HOA.pdf)
    # we must account as well for the Condon-Shortley phase, which is already included in coef
    sn3d_factor = pow(-1,m) *np.sqrt ( (2-delta_kronecker(0,m)) * ( float(factorial(l-m)) / float(factorial(l+m)) ) )
    coef = coef / get_spherical_harmonic_normalization_coef(ambisonics_order,ambisonics_degree) * sn3d_factor

    # if we would want to include N3D support at some moment, uncomment this lines
    n3d_factor = np.sqrt((2*l)+1)
    coef = coef * n3d_factor

    # TODO: COMPATIBILITY WITH MATLAB CODE...
    # this is approx 1/(2*sqrt(3)), but still don't know where it comes from...
    coef = coef * 0.28209479177387814

    return coef


def get_ambisonics_coefs(azimuth,elevation,order):

    coefs = np.zeros(np.power(order+1,2))
    coef_index = 0
    for l in range(order+1):
        for m in range(-l,l+1): # this is ACN channel ordering
            coefs[coef_index] = get_spherical_harmonic(azimuth, elevation, l, m)
            coef_index += 1
    return coefs



