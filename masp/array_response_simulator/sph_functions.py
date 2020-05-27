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
from scipy.special import spherical_jn, spherical_yn
from masp.validate_data_types import _validate_int, _validate_ndarray_1D

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Spherical Hankel functions

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
    return spherical_jn(n, x) + 1j * spherical_yn(n, x)

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
    return spherical_jn(n, x) - 1j * spherical_yn(n, x)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# First derivative of spherical Hankel functions

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
    return spherical_jn(n, x, derivative=True) + 1j * spherical_yn(n, x, derivative=True)

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
    return spherical_jn(n, x, derivative=True) - 1j * spherical_yn(n, x, derivative=True)
