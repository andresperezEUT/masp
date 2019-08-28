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
#   @file   validate_data_types.py
#   @author Andrés Pérez-López
#   @date   31/07/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import numpy as np
from masp.shoebox_room_sim.echogram import Echogram

def _validate_boolean(name, boolean):

    if not isinstance(boolean, bool):
        raise TypeError(name + ' must be an instance of bool')

def _validate_int(name, number, positive=False, limit=None, parity=None):

    if not isinstance(number, int):
        raise TypeError(name + ' must be an instance of int')
    if positive:
        if number < 0:
            raise ValueError(name + ' must be positive')
    if limit is not None:
        if number > limit:
            raise ValueError(name + ' must be smaller than ' + str(limit))
    if parity is not None:
        if parity == 'even':
            if number%2 != 0:
                raise ValueError(name + ' must be even')
        elif parity == 'odd':
            if number%2 == 0:
                raise ValueError(name + ' must be odd')
        else:
            raise TypeError(name + ': unknown parity value: ' + parity)

def _validate_float(name, number, positive=False, limit=None):

    if not isinstance(number, float):
        raise TypeError(name + ' must be an instance of float')
    if positive:
        if number < 0:
            raise ValueError(name + ' must be positive')
    if limit is not None:
        if number > limit:
            raise ValueError(name + ' must be smaller than ' + str(limit))

def _validate_number(name, number, norm=False, positive=False, limit=None):

    if isinstance(number, (int, float)):
        if norm:
            if not 0 <= number <= 1:
                raise ValueError(name + ' must be in the interval [0,1]')
        if positive:
            if 0 > number:
                raise ValueError(name + ' must be positive')
        if limit is not None:
            if number < limit[0]:
                raise ValueError(name + ' must be greater or equal than ' + str(limit[0]))
            if number > limit[1]:
                raise ValueError(name + ' must be smaller or equal than ' + str(limit[1]))
    elif isinstance(number, np.ndarray):
        if norm:
            if np.any(number < 0) or np.any(number > 1):
                raise ValueError(name + ' must be in the interval [0,1]')
        if positive:
            if np.any(number < 0):
                raise ValueError(name + ' must be positive')
        if limit is not None:
            if np.any(number < limit[0]):
                raise ValueError(name + ' must be greater or equal than ' + str(limit[0]))
            if np.any(number > limit[1]):
                raise ValueError(name + ' must be smaller or equal than ' + str(limit[1]))
    else:
        raise TypeError(name + ' must be an integer, float or 1-D ndarray')

def _validate_ndarray_1D(name, ndarray, size=None, norm=False, positive=False, limit=None, dtype=None):

    if not isinstance(ndarray, np.ndarray):
        raise TypeError(name + ' must be an instance of ndarray')
    if ndarray.ndim != 1:
        raise ValueError(name + ' must be 1D')
    if size is not None:
        if ndarray.size != size:
            raise ValueError(name + ' must have size ' + str(size))
    if norm:
        if np.any(ndarray < 0) or np.any(ndarray > 1):
            raise ValueError(name + ' must be in the interval [0,1]')
    if positive:
        if np.any(ndarray < 0):
            raise ValueError(name + ' must be positive')
    if limit is not None:
        if np.any(ndarray < limit[0]):
            raise ValueError(name + ' must be greater or equal than ' + str(limit[0]))
        if np.any(ndarray > limit[1]):
            raise ValueError(name + ' must be smaller or equal than ' + str(limit[1]))
    if dtype is not None:
        if ndarray.dtype != dtype:
            raise ValueError(name + ' dtype must be ' + str(dtype))

def _validate_ndarray_2D(name, ndarray, shape0=None, shape1=None, norm=False, positive=False, dtype=None):

    if not isinstance(ndarray, np.ndarray):
        raise TypeError(name + ' must be an instance of ndarray')
    if ndarray.ndim != 2:
        raise ValueError(name + ' must be 2D')
    if shape0 is not None:
        if ndarray.shape[0] != shape0:
            raise ValueError(name + ' must have dimension 0=' + str(shape0))
    if shape1 is not None:
        if ndarray.shape[1] != shape1:
            raise ValueError(name + ' must have dimension 1=' + str(shape1))
    if norm:
        if np.any(ndarray < 0) or np.any(ndarray > 1):
            raise ValueError(name + ' must be in the interval [0,1]')
    if positive:
        if np.any(ndarray < 0):
            raise ValueError(name + ' must be positive')
    if dtype is not None:
        if ndarray.dtype != dtype:
            raise ValueError(name + ' dtype must be ' + str(dtype))

def _validate_ndarray_3D(name, ndarray, shape0=None, shape1=None, shape2=None, norm=False, positive=False, dtype=None):

    if not isinstance(ndarray, np.ndarray):
        raise TypeError(name + ' must be an instance of ndarray')
    if ndarray.ndim != 3:
        raise ValueError(name + ' must be 3D')
    if shape0 is not None:
        if ndarray.shape[0] != shape0:
            raise ValueError(name + ' must have dimension 0=' + str(shape0))
    if shape1 is not None:
        if ndarray.shape[1] != shape1:
            raise ValueError(name + ' must have dimension 1=' + str(shape1))
    if shape2 is not None:
        if ndarray.shape[2] != shape2:
            raise ValueError(name + ' must have dimension 2=' + str(shape1))
    if norm:
        if np.any(ndarray < 0) or np.any(ndarray > 1):
            raise ValueError(name + ' must be in the interval [0,1]')
    if positive:
        if np.any(ndarray < 0):
            raise ValueError(name + ' must be positive')
    if dtype is not None:
        if ndarray.dtype != dtype:
            raise ValueError(name + ' dtype must be ' + str(dtype))

def _validate_ndarray_4D(name, ndarray, shape0=None, shape1=None, shape2=None, shape3=None, norm=False, positive=False, dtype=None):

    if not isinstance(ndarray, np.ndarray):
        raise TypeError(name + ' must be an instance of ndarray')
    if ndarray.ndim != 4:
        raise ValueError(name + ' must be 4D')
    if shape0 is not None:
        if ndarray.shape[0] != shape0:
            raise ValueError(name + ' must have dimension 0=' + str(shape0))
    if shape1 is not None:
        if ndarray.shape[1] != shape1:
            raise ValueError(name + ' must have dimension 1=' + str(shape1))
    if shape2 is not None:
        if ndarray.shape[2] != shape2:
            raise ValueError(name + ' must have dimension 2=' + str(shape1))
    if shape3 is not None:
        if ndarray.shape[3] != shape2:
            raise ValueError(name + ' must have dimension 3=' + str(shape1))
    if norm:
        if np.any(ndarray < 0) or np.any(ndarray > 1):
            raise ValueError(name + ' must be in the interval [0,1]')
    if positive:
        if np.any(ndarray < 0):
            raise ValueError(name + ' must be positive')
    if dtype is not None:
        if ndarray.dtype != dtype:
            raise ValueError(name + ' dtype must be ' + str(dtype))

def _validate_string(name, string, choices=None):

    if not isinstance(string, str):
        raise TypeError(name + ' must be an instance of str')
    if choices is not None and string not in choices:
        raise ValueError(name + ' must be one of the following: ' + str(choices))


#
# def _validate_room(room):
#     if not isinstance(room, np.ndarray):
#         raise TypeError('room must be an instance of ndarray')
#     elif room.ndim != 1 or room.size != 3:
#         raise ValueError('room must have dimension (3)')
#     elif np.any(room < 0):
#         raise ValueError('room coordinates must be positive')
#
# def _validate_rt60_target(rt60_target):
#     if not isinstance(rt60_target, np.ndarray):
#         raise TypeError('rt60_target must be an instance of ndarray')
#     elif rt60_target.ndim != 1:
#         raise ValueError('rt60_target must have dimension (nBands)')
#     elif np.any(rt60_target < 0):
#         raise ValueError('rt60_target values must be positive')
#
# def _validate_abs_wall_ratios(abs_wall_ratios, force_normalized=False):
#     if not isinstance(abs_wall_ratios, np.ndarray):
#         raise TypeError('abs_wall_ratios must be an instance of ndarray, or None')
#     elif abs_wall_ratios.ndim != 1 or abs_wall_ratios.size != 6:
#         raise ValueError('abs_wall_ratios must have dimension (6)')
#     if force_normalized:
#         if np.any(abs_wall_ratios < 0) or np.any(abs_wall_ratios > 1):
#             raise ValueError('abs_wall_ratios must be in the interval [0,1]')
#
# def _validate_alpha(alpha):
#
#     if isinstance(alpha, (int, float)):
#         if alpha < 0 or alpha > 1:
#             raise ValueError('alpha must be in the interval [0,1]')
#     elif isinstance(alpha, np.ndarray):
#         if alpha.shape != (1,):
#             raise TypeError('if alpha is ndarray, it must be unidimensional')
#         elif np.any(alpha < 0) or np.any(alpha > 1):
#             raise ValueError('alpha must be in the interval [0,1]')
#     else:
#         raise TypeError('alpha must be an integer, float or 1-D ndarray')

def _validate_echogram(echogram):
    from masp.utils import C

    if not isinstance(echogram, Echogram):
        raise TypeError('echogram must be an instance of Echogram')

    _validate_ndarray_2D('echogram.value', echogram.value)
    _validate_ndarray_1D('echogram.time', echogram.time, positive=True)
    _validate_ndarray_2D('echogram.order', echogram.order, shape1=C, dtype=int)
    _validate_ndarray_2D('echogram.coords', echogram.coords, shape1=C)

    if echogram.value.shape[0] != echogram.time.shape[0] != echogram.order.shape[0] != echogram.coords.shape[0]:
        raise ValueError('echogram shape mismatch')

def _validate_echogram_array(ndarray, shape0=None, shape1=None, shape2=None):
    """
    specific case of 2D/3D ndarray with dtype=Echogram
    """
    if not isinstance(ndarray, np.ndarray):
        raise TypeError('Echogram array must be an instance of ndarray')
    if ndarray.ndim not in (2, 3):
        raise ValueError('Echogram array must be 2D or 3D')
    if shape0 is not None:
        if ndarray.shape[0] != shape0:
            raise ValueError('Echogram array must have dimension 0=' + str(shape0))
    if shape1 is not None:
        if ndarray.shape[1] != shape1:
            raise ValueError('Echogram array must have dimension 1=' + str(shape1))
    if shape2 is not None:
        if ndarray.shape[2] != shape1:
            raise ValueError('Echogram array must have dimension 2=' + str(shape2))
    if ndarray.dtype != Echogram:
        raise ValueError('Echogram array dtype must be Echogram')
    for idx in np.ndindex(ndarray.shape):
        _validate_echogram(ndarray[idx])

# def _validate_alpha_walls_per_band(alpha):
#     if not isinstance(alpha, np.ndarray):
#         raise TypeError('alpha must be an instance of ndarray')
#     elif alpha.ndim != 2 or alpha.shape[1] != 6:
#         raise ValueError('alpha must have dimension (nBands, 6)')
#     elif np.any(alpha < 0) or np.any(alpha > 1):
#         raise ValueError('alpha values must be in the interval [0,1]')
#
#
# def _validate_time_limits(limits, nBands):
#     if not isinstance(limits, np.ndarray):
#         raise TypeError('limits must be an instance of ndarray, or None')
#     elif limits.ndim != 1 or limits.size != nBands:
#         raise TypeError('limits dimension must be (nBands,)')