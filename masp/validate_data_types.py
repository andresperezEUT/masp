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
from masp.shoebox_room_sim.echogram import Echogram, QuantisedEchogram

def _validate_boolean(name, boolean):

    if not isinstance(boolean, bool):
        raise TypeError(name + ' must be an instance of bool')

def _validate_int(name, number, positive=False, limit=None, parity=None):

    if not isinstance(number, int) or isinstance(number, bool):
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

    if isinstance(number, (int, float)) and not isinstance(number, bool):
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

def _validate_list(name, arg, size=None):

    if not isinstance(arg, list):
        raise TypeError(name + ' must be an instance of list')
    if size is not None:
        if len(arg) != size:
            raise ValueError(name + ' must have size ' + str(size))

def _validate_ndarray(name, ndarray):

    if not isinstance(ndarray, np.ndarray):
        raise TypeError(name + ' must be an instance of ndarray')

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
            raise TypeError(name + ' dtype must be ' + str(dtype))

def _validate_ndarray_2D(name, ndarray, shape0=None, shape1=None, norm=False, positive=False, limit=None, dtype=None):

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
    if limit is not None:
        if np.any(ndarray < limit[0]):
            raise ValueError(name + ' must be greater or equal than ' + str(limit[0]))
        if np.any(ndarray > limit[1]):
            raise ValueError(name + ' must be smaller or equal than ' + str(limit[1]))
    if dtype is not None:
        if ndarray.dtype != dtype:
            raise TypeError(name + ' dtype must be ' + str(dtype))

def _validate_ndarray_3D(name, ndarray, shape0=None, shape1=None, shape2=None, norm=False, positive=False, limit=None, dtype=None):

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
    if limit is not None:
        if np.any(ndarray < limit[0]):
            raise ValueError(name + ' must be greater or equal than ' + str(limit[0]))
        if np.any(ndarray > limit[1]):
            raise ValueError(name + ' must be smaller or equal than ' + str(limit[1]))
    if dtype is not None:
        if ndarray.dtype != dtype:
            raise TypeError(name + ' dtype must be ' + str(dtype))

def _validate_ndarray_4D(name, ndarray, shape0=None, shape1=None, shape2=None, shape3=None, norm=False, positive=False, limit=None, dtype=None):

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
    if limit is not None:
        if np.any(ndarray < limit[0]):
            raise ValueError(name + ' must be greater or equal than ' + str(limit[0]))
        if np.any(ndarray > limit[1]):
            raise ValueError(name + ' must be smaller or equal than ' + str(limit[1]))
    if dtype is not None:
        if ndarray.dtype != dtype:
            raise TypeError(name + ' dtype must be ' + str(dtype))

def _validate_string(name, string, choices=None):

    if not isinstance(string, str):
        raise TypeError(name + ' must be an instance of str')
    if choices is not None and string not in choices:
        raise ValueError(name + ' must be one of the following: ' + str(choices))


def _validate_echogram(echogram):
    from masp.utils import C

    if not isinstance(echogram, Echogram):
        raise TypeError('echogram must be an instance of Echogram')

    _validate_ndarray_2D('echogram.value', echogram.value)
    _validate_ndarray_1D('echogram.time', echogram.time, positive=True)
    _validate_ndarray_2D('echogram.order', echogram.order, shape1=C, dtype=int)
    _validate_ndarray_2D('echogram.coords', echogram.coords, shape1=C)

    shapes = [echogram.value.shape[0], echogram.time.shape[0], echogram.order.shape[0],  echogram.coords.shape[0]]
    if not all(s == shapes[0] for s in shapes):
        raise ValueError('echogram shape mismatch')


def _validate_quantised_echogram(qechogram):
    from masp.utils import C

    if not isinstance(qechogram, QuantisedEchogram):
        raise TypeError('quantised echogram must be an instance of QuantisedEchogram')

    _validate_ndarray_2D('qechogram.value', qechogram.value)
    _validate_ndarray_1D('qechogram.time', qechogram.time, positive=True)
    _validate_boolean('qechogram.isActive', qechogram.isActive)

    if qechogram.value.shape[0] != qechogram.time.shape[0]:
        raise ValueError('quantised echogram shape mismatch')

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
    if ndarray.dtype != Echogram or not isinstance(ndarray.flatten()[0], Echogram):
        raise TypeError('Echogram array dtype must be Echogram')
    for idx in np.ndindex(ndarray.shape):
        _validate_echogram(ndarray[idx])


def _validate_quantised_echogram_array(ndarray, size=None):
    """
    specific case of 1D ndarray with dtype=QuantisedEchogram
    """
    if not isinstance(ndarray, np.ndarray):
        raise TypeError('Quantised Echogram array must be an instance of ndarray')
    if ndarray.ndim != 1:
        raise ValueError('Quantised Echogram array must be 1D')
    if size is not None:
        if ndarray.size != size:
            raise ValueError('Quantised Echogram array must have size=' + str(size))
    if ndarray.dtype != QuantisedEchogram  or not isinstance(ndarray.flatten()[0], QuantisedEchogram):
        raise TypeError('Quantised Echogram array dtype must be QuantisedEchogram')
    for idx in np.ndindex(ndarray.shape):
        _validate_quantised_echogram(ndarray[idx])
