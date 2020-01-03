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
import csv
import numpy as np
from scipy.special import sph_harm
from masp.validate_data_types import _validate_int, _validate_ndarray_2D, _validate_string, _validate_ndarray_1D, \
    _validate_ndarray

c = 343.
C = 3

def get_capsule_positions(mic_array_name):
    """
    Retrieve the geometry of a selected set of microphone arrays.

    Parameters
    ----------
    mic_array_name : str
      One of: 'eigenmike', 'ambeo'

    Returns
    -------
    array_sigs : ndarray
        Capsule positions, in spherical coordinates (radians). Dimension = (nMic, C)

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    """

    _validate_string('mic_array_name', mic_array_name, choices=['eigenmike', 'ambeo'])
    capsule_positions = None

    if mic_array_name is 'eigenmike':
        mic_dirs_deg = np.array([[0, 32, 0, 328, 0, 45, 69, 45, 0, 315, 291, 315, 91, 90, 90, 89, 180, 212, 180, 148,
                                  180, 225, 249, 225, 180, 135, 111, 135, 269, 270, 270, 271],
                                 [21, 0, -21, 0, 58, 35, 0, -35, -58, -35, 0, 35, 69, 32, -31, -69, 21, 0, -21, 0, 58,
                                  35, 0, -35, -58, -35, 0, 35, 69, 32, -32, -69]])
        mic_dirs_rad = mic_dirs_deg * np.pi / 180.
        r = 0.042
        mic_dirs_rad = np.row_stack((mic_dirs_rad, r*np.ones(np.shape(mic_dirs_rad)[1])))

        capsule_positions = mic_dirs_rad.T

    elif mic_array_name is 'ambeo':
        r = 0.015
        capsule_positions = [[np.pi / 4, np.arcsin(1. / np.sqrt(3)), r],           # FLU
                             [7 * np.pi / 4, -1 * np.arcsin(1. / np.sqrt(3)), r],  # FRD
                             [3 * np.pi / 4, -1 * np.arcsin(1. / np.sqrt(3)), r],  # BLD
                             [5 * np.pi / 4, np.arcsin(1. / np.sqrt(3)), r]]       # BRU

    return capsule_positions



def cart2sph(cart):
    """
    Cartesian to spherical coordinates transformation, in matrix form.

    Parameters
    ----------
    cart : ndarray
      Cartesian coordinates. Dimension = (nCoords, C)

    Returns
    -------
    sph : ndarray
        Spherical coordinates, in radians, aed.  Dimension = (nCoords, C)

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    As a dimensionality exception, in case the input matrix is 1D (just one point),
    the output matrix will be as well 1D.
    """

    arg = cart.copy()
    _validate_ndarray('cart', cart)
    if cart.ndim == 1:
        _validate_ndarray_1D('cart', cart, size=C)
        cart = cart[np.newaxis, :]
    elif cart.ndim == 2:
        _validate_ndarray_2D('cart', cart, shape1=C)
    else:
        raise ValueError('cart must be either 1D or 2D array')

    sph = np.empty(cart.shape)
    hypotxy = np.hypot(cart[:, 0], cart[:, 1])
    sph[:, 2] = np.hypot(hypotxy, cart[:, 2])
    sph[:, 1] = np.arctan2(cart[:, 2], hypotxy)
    sph[:, 0] = np.arctan2(cart[:, 1], cart[:, 0])
    if arg.ndim == 1:
        sph = sph.squeeze()
    return sph


def sph2cart(sph):
    """
    Spherical to cartesian coordinates transformation, in matrix form.

    Parameters
    ----------
    sph : ndarray
        Spherical coordinates, in radians, aed.  Dimension = (nCoords, C)


    Returns
    -------
    sph : ndarray
        Cartesian coordinates. Dimension = (nCoords, C)

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    As a dimensionality exception, in case the input matrix is 1D (just one point),
    the output matrix will be as well 1D.
    """

    arg = sph.copy()
    _validate_ndarray('sph', sph)
    if sph.ndim == 1:
        _validate_ndarray_1D('sph', sph, size=C)
        sph = sph[np.newaxis, :]
    elif sph.ndim == 2:
        _validate_ndarray_2D('sph', sph, shape1=C)
    else:
        raise ValueError('sph must be either 1D or 2D array')

    cart = np.empty(sph.shape)
    cart[:, 2] = sph[:, 2] * np.sin( sph[:, 1])
    rcoselev = sph[:, 2] * np.cos( sph[:, 1])
    cart[:, 0] = rcoselev * np.cos( sph[:, 0])
    cart[:, 1] = rcoselev * np.sin( sph[:, 0])
    if arg.ndim == 1:
        cart = cart.squeeze()
    return cart


def elev2incl(dirs):
    """
    Spherical coordinates: elevation to inclination reference system

    Parameters
    ----------
    dirs : ndarray
        Spherical coordinates, in radians, aed.  Dimension = (nCoords, {2,3})


    Returns
    -------
    incl : ndarray
       Transformed coordinates. Dimension = (nCoords, {2,3})

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    The input matrix might have dimension 1 = 2 ([azimuth, elevation]),
    or dimension 1 = 3 ([azimuth, elevation, distance]).
    The output matrix will propagate input dimensionality.
    """

    _validate_ndarray_2D('dirs',dirs)
    if dirs.shape[1] == 2:
        incl = np.column_stack((dirs[:, 0], np.pi / 2 - dirs[:, 1]))
    elif dirs.shape[1] == 3:
        incl = np.column_stack((dirs[:, 0], np.pi / 2 - dirs[:, 1], dirs[:, 2]))
    else:
        raise ValueError('dirs must have dimension 1={2,3}')
    return incl


def incl2elev(dirs):
    """
    Spherical coordinates: inclination to elevation reference system

    Parameters
    ----------
    dirs : ndarray
        Spherical coordinates, in radians, aid.  Dimension = (nCoords, {2,3})

    Returns
    -------
    elev : ndarray
       Transformed coordinates. Dimension = (nCoords, {2,3})

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    The input matrix might have dimension 1 = 2 ([azimuth, inclination]),
    or dimension 1 = 3 ([azimuth, inclination, distance]).
    The output matrix will propagate input dimensionality.
    """
    return elev2incl(dirs)




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


def isLambda(v):
    """
    Determine if a given argument is a lambda expression.

    Parameters
    ----------
    v : arg
        Argument to test.

    Returns
    -------
    isLambda : boolean
        Result.
    """

    LAMBDA = lambda:0
    return isinstance(v, type(LAMBDA)) and v.__name__ == LAMBDA.__name__


def load_sph_grid(file_path):
    """
    todo: implementation_ at the moment only 2D matrix
    :param file_path:
    :return:
    """

    with open(file_path) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=' ')
        rows = list(csv_reader)
        dim1 = int(rows[0][0])
        dim0 = int(rows[1][0])
        sph_grid = np.asarray(rows[2:], dtype='float')

    assert np.shape(sph_grid) == (dim0, dim1)
    return sph_grid



#############################################################################################
#
# Spherical Harmonic Transform library.
# TODO: separate package?

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
        Spherical harmonics . Dimension = (nDirs, nHarm).

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
    _validate_ndarray_2D('dirs', dirs, shape1=C-1)
    _validate_string('basisType', basisType, choices=['complex', 'real'])

    nDirs = dirs.shape[0]
    nHarm = np.power(N+1, 2)
    Y_N = np.zeros((nDirs, nHarm))

    def delta_kronecker(q1, q2):
        return 1 if q1 == q2 else 0

    # TODO
    # it looks like the output of shs is N3d (max 1, sqrt(3)!3)
    # so it needs to be scaled as * np.sqrt(4*np.pi) * [1, 1./np.sqrt(3), 1./np.sqrt(3), 1./np.sqrt(3)]

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
        # TODO
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


def check_cond_number_sht(N, dirs, basisType, W=None):
    """
    Computes the condition number for a least-squares SHT.

    Parameters
    ----------
    N : int
        Maximum order to be tested for the given set of points.
   dirs : ndarray
        Evaluation directions. Dimension = (nDirs, 2).
        Directions are expected in radians, expressed in pairs [azimuth, inclination].
    basisType : str
        Type of spherical harmonics. Either 'complex' or 'real'.
    W : ndarray, optional.
        Weights for each measurement point to condition the inversion,
        in case of a weighted least-squares. Dimension = (nDirs)

    Returns
    -------
    cond_N : ndarray
        Condition number for each order. Dimension = (N+1)

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    Inclination is defined as the angle from zenith: inclination = pi/2-elevation

    TODO: implement complex basis?
    TODO: implement W
    """

    _validate_int('N', N, positive=True)
    _validate_ndarray_2D('dirs', dirs, shape1=C-1)
    _validate_string('basisType', basisType, choices=['complex', 'real'])
    if W is not None:
        _validate_ndarray_1D('W', W, size=dirs.shape[0])

    # Compute the harmonic coefficients
    Y_N = get_sh(N, dirs, basisType)

    # Compute condition number for progressively increasing order up to N
    cond_N = np.zeros(N + 1)
    for n in range(N+1):
        Y_n = Y_N[:, :np.power(n + 1, 2)]
        if W is None:
            YY_n = np.dot(Y_n.T, Y_n)
        else:
            # YY_n = Y_n.T * np.diag(W) * Y_n
            raise NotImplementedError
        cond_N[n] = np.linalg.cond(YY_n)

    return cond_N


def replicate_per_order(x):
    """
    Replicate l^th element 2*l+1 times across dimension

    Parameters
    ----------
    x : ndarray
        Array to replicate. Dimension = (l)

    Returns
    -------
    x_rep : ndarray
        Replicated array. Dimension = ( (l+1)^2 )

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    Replicates multidimensional array across dimension dim, so that the
    l^th element of dim is replicated 2*l+1 times. that effectively has the
    effect that the dimension grows from L to L^2 elements. This can be useful
    in some spherical harmonic operations.

    TODO: FOR THE MOMENT JUST 1D
    TODO: optimize
    """

    _validate_ndarray_1D('x', x)

    order = np.size(x) - 1
    n_sh = np.power(order+1, 2)
    x_rep = np.zeros(n_sh, dtype=x.dtype)
    sh_idx = 0

    for m in range(order+1):
        n_sh_order = 2 * m + 1  # number of spherical harmonics at the given order m
        for n in range(n_sh_order):
            x_rep[sh_idx] = x[m]
            sh_idx += 1

    return x_rep
