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
#   @file   sph_functions.py
#   @author Andrés Pérez-López
#   @date   22/08/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import numpy as np
from scipy.special import jv, yv
from masp.validate_data_types import _validate_int, _validate_ndarray_1D

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Spherical functions


def sph_besselj(n, x):
    """
    Spherical bessel function of the first kind.

    Parameters
    ----------
    n : int
        Function order.
    x: ndarray
        Points where to evaluate the function. Dimension = (l)

    Returns
    -------
    f : ndarray
        Function result. Dimension = (l)

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    """
    _validate_int('n', n)
    _validate_ndarray_1D('x', x)
    idx_zero, = np.where(x == 0)
    j = np.sqrt(np.pi / (2 * x)) * jv(n + 0.5, x)
    j[idx_zero] = 1. if n == 0 else 0.
    return j

def sph_bessely(n, x):
    """
    Spherical bessel function of the second kind.

    Parameters
    ----------
    n : int
        Function order.
    x: ndarray
        Points where to evaluate the function. Dimension = (l)

    Returns
    -------
    f : ndarray
        Function result. Dimension = (l)

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    """
    _validate_int('n', n)
    _validate_ndarray_1D('x', x)
    return np.sqrt(np.pi / (2 * x)) * yv(n + 0.5, x)

def sph_function():
    raise NotImplementedError

def sph_hankel1(n, x):
    """
    Spherical hankel function of the first kind.

    Parameters
    ----------
    n : int
        Function order.
    x: ndarray
        Points where to evaluate the function. Dimension = (l)

    Returns
    -------
    f : ndarray
        Function result. Dimension = (l)

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    """
    _validate_int('n', n)
    _validate_ndarray_1D('x', x)
    return sph_besselj(n, x) + 1j * sph_bessely(n, x)

def sph_hankel2(n, x):
    """
    Spherical hankel function of the second kind.

    Parameters
    ----------
    n : int
        Function order.
    x: ndarray
        Points where to evaluate the function. Dimension = (l)

    Returns
    -------
    f : ndarray
        Function result. Dimension = (l)

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    """
    _validate_int('n', n)
    _validate_ndarray_1D('x', x)
    return sph_besselj(n, x) - 1j * sph_bessely(n, x)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# First derivative of spherical functions

def dsph_besselj(n, x):
    """
    Spherical bessel function derivative of the first kind.

    Parameters
    ----------
    n : int
        Function order.
    x: ndarray
        Points where to evaluate the function. Dimension = (l)

    Returns
    -------
    f : ndarray
        Function result. Dimension = (l)

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    """
    _validate_int('n', n)
    _validate_ndarray_1D('x', x)
    return 1. / (2 * n + 1) * (n * sph_besselj(n - 1, x) - (n + 1) * sph_besselj(n + 1, x))

def dsph_bessely(n, x):
    """
    Spherical bessel function derivative of the second kind.

    Parameters
    ----------
    n : int
        Function order.
    x: ndarray
        Points where to evaluate the function. Dimension = (l)

    Returns
    -------
    f : ndarray
        Function result. Dimension = (l)

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    """
    _validate_int('n', n)
    _validate_ndarray_1D('x', x)
    return 1. / (2 * n+1) * (n*sph_bessely(n-1, x) - (n+1)*sph_bessely(n+1, x))

def dsph_hankel1(n, x):
    """
    Spherical hankel function derivative of the first kind.

    Parameters
    ----------
    n : int
        Function order.
    x: ndarray
        Points where to evaluate the function. Dimension = (l)

    Returns
    -------
    f : ndarray
        Function result. Dimension = (l)

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    """
    _validate_int('n', n)
    _validate_ndarray_1D('x', x)
    return dsph_besselj(n, x) + 1j * dsph_bessely(n, x)

def dsph_hankel2(n, x):
    """
    Spherical hankel function derivative of the second kind.

    Parameters
    ----------
    n : int
        Function order.
    x: ndarray
        Points where to evaluate the function. Dimension = (l)

    Returns
    -------
    f : ndarray
        Function result. Dimension = (l)

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    """
    _validate_int('n', n)
    _validate_ndarray_1D('x', x)
    return dsph_besselj(n, x) - 1j * dsph_bessely(n, x)
