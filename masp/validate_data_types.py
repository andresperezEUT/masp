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

def _validate_int(name, number, positive=False, limit=None):

    if not isinstance(number, int):
        raise TypeError(name + ' must be an instance of int')
    elif positive:
        if number < 0:
            raise ValueError(name + ' must be positive')
    elif limit is not None:
        if number > limit:
            raise ValueError(name + ' must be smaller than ' + str(limit))

def _validate_number(name, number, norm=False, positive=False, limit=None):

    if isinstance(number, (int, float)):
        if norm:
            if not 0 <= number <= 1:
                raise ValueError(name + ' must be in the interval [0,1]')
        elif positive:
            if 0 > number:
                raise ValueError(name + ' must be positive')
        elif limit is not None:
            if number > limit:
                raise ValueError(name + ' must be smaller than ' + str(limit))
    elif isinstance(number, np.ndarray):
        if norm:
            if np.any(number < 0) or np.any(number > 1):
                raise ValueError(name + ' must be in the interval [0,1]')
        elif positive:
            if np.any(number < 0):
                raise ValueError(name + ' must be positive')
        elif limit is not None:
            if np.any(number > limit):
                raise ValueError(name + ' must be smaller than ' + str(limit))
    else:
        raise TypeError(name + ' must be an integer, float or 1-D ndarray')

def _validate_ndarray_1D(name, ndarray, size=None, norm=False, positive=False, limit=None, dtype=None):

    if not isinstance(ndarray, np.ndarray):
        raise TypeError(name + ' must be an instance of ndarray')
    elif ndarray.ndim != 1:
        raise ValueError(name + ' must be 1D')
    elif size is not None:
        if ndarray.size != size:
            raise ValueError(name + ' must have size ' + str(size))
    elif norm:
        if np.any(ndarray < 0) or np.any(ndarray > 1):
            raise ValueError(name + ' must be in the interval [0,1]')
    elif positive:
        if np.any(ndarray < 0):
            raise ValueError(name + ' must be positive')
    elif limit is not None:
        if np.any(ndarray > limit):
            raise ValueError(name + ' must be smaller than ' + str(limit))
    elif dtype is not None:
        if ndarray.dtype != dtype:
            raise ValueError(name + ' dtype must be ' + dtype)

def _validate_ndarray_2D(name, ndarray, shape0=None, shape1=None, norm=False, positive=False, dtype=None):

    if not isinstance(ndarray, np.ndarray):
        raise TypeError(name + ' must be an instance of ndarray')
    elif ndarray.ndim != 2:
        raise ValueError(name + ' must be 2D')

    elif shape0 is not None:
        if ndarray.shape[0] != shape0:
            raise ValueError(name + ' must have dimension 0=' + str(shape0))
    elif shape1 is not None:
        if ndarray.shape[1] != shape1:
            raise ValueError(name + ' must have dimension 1=' + str(shape1))
    elif norm:
        if np.any(ndarray < 0) or np.any(ndarray > 1):
            raise ValueError(name + ' must be in the interval [0,1]')
    elif positive:
        if np.any(ndarray < 0):
            raise ValueError(name + ' must be positive')
    elif dtype is not None:
        if ndarray.dtype != dtype:
            raise ValueError(name + ' dtype must be ' + dtype)

def _validate_string(name, string, choices=None):

    if not isinstance(string, str):
        raise TypeError(name + ' must be an instance of str')

    elif choices is not None and string not in choices:
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

    if not isinstance(echogram, Echogram):
        raise TypeError('echogram must be an instance of Echogram')

    _validate_ndarray_1D('echogram.value', echogram.value, positive=True)
    # elif not isinstance(echogram.value, np.ndarray):
    #     raise TypeError('echogram.value must be a ndarray')
    # elif echogram.value.ndim != 1:
    #     raise ValueError('echogram.value dimension must be (t,)')

    _validate_ndarray_1D('echogram.time', echogram.time, positive=True)
    # elif not isinstance(echogram.time, np.ndarray):
    #     raise TypeError('echogram.time must be a ndarray')
    # elif echogram.time.ndim != 1:
    #     raise ValueError('echogram.time dimension must be (t,)')

    _validate_ndarray_2D('echogram.order', echogram.order, shape1=3)
    # order values should be integers
    if not np.all(echogram.order == np.floor(echogram.order)):
        raise ValueError('echogram.order should contain integers')

    # elif not isinstance(echogram.order, np.ndarray):
    #     raise TypeError('echogram.order must be a ndarray')
    # elif echogram.order.ndim != 2 or echogram.order.shape[1] !=3:
    #     raise ValueError('echogram.order dimension must be (t, 3)')

    _validate_ndarray_2D('echogram.coords', echogram.coords, shape1=3)
    # elif not isinstance(echogram.coords, np.ndarray):
    #     raise TypeError('echogram.coords must be a ndarray')
    # elif echogram.coords.ndim != 2 or echogram.coords.shape[1] !=3:
    #     raise ValueError('echogram.coords dimension must be (t, 3)')

    if echogram.value.shape[0] != echogram.time.shape[0] != echogram.order.shape[0] != echogram.coords.shape[0]:
        raise ValueError('echogram shape mismatch')

def _validate_echogram_array(ndarray, shape0=None, shape1=None):
    """
    specific case of 2D ndarray with dtype=Echogram
    """
    if not isinstance(ndarray, np.ndarray):
        raise TypeError('Echogram array must be an instance of ndarray')
    elif ndarray.ndim != 2:
        raise ValueError('Echogram array must be 2D')
    elif shape0 is not None:
        if ndarray.shape[0] != shape0:
            raise ValueError('Echogram array must have dimension 0=' + str(shape0))
    elif shape1 is not None:
        if ndarray.shape[1] != shape1:
            raise ValueError('Echogram array must have dimension 1=' + str(shape1))
    elif ndarray.dtype != Echogram:
        raise ValueError('Echogram array dtype must be Echogram')
    for i in range(ndarray.shape[0]):
        for j in range(ndarray.shape[1]):
            _validate_echogram(ndarray[i,j])

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