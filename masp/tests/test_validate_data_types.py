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
#   @file   test_validate_data_types.py
#   @author Andrés Pérez-López
#   @date   25/09/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import pytest
import numpy as np

from asma import Echogram, QuantisedEchogram, C
from asma.validate_data_types import _validate_boolean
from asma.validate_data_types import _validate_int
from asma.validate_data_types import _validate_float
from asma.validate_data_types import _validate_number
from asma.validate_data_types import _validate_list
from asma.validate_data_types import _validate_ndarray
from asma.validate_data_types import _validate_ndarray_1D
from asma.validate_data_types import _validate_ndarray_2D
from asma.validate_data_types import _validate_ndarray_3D
from asma.validate_data_types import _validate_ndarray_4D
from asma.validate_data_types import _validate_string
from asma.validate_data_types import _validate_echogram
from asma.validate_data_types import _validate_quantised_echogram
from asma.validate_data_types import _validate_echogram_array
from asma.validate_data_types import _validate_quantised_echogram_array


def test_validate_boolean():

    # TypeError: not a boolean
    wrong_values = ['True', 1, 2.3, 1e4, 3j, [True], None, np.nan, np.inf, np.asarray([0.5])]
    for wv in wrong_values:
        with pytest.raises(TypeError, match='must be an instance of bool'):
            _validate_boolean('vw', wv)


def test_validate_int():

    # TypeError: not an integer
    wrong_values = ['1', True, 2.3, 1e4, 3j, [1], None, np.nan, np.inf, np.asarray([0.5])]
    for wv in wrong_values:
        with pytest.raises(TypeError, match='must be an instance of int'):
            _validate_int('vw', wv)

    # ValueError: not positive
    wrong_values = [-1, -3]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must be positive'):
            _validate_int('vw', wv, positive=True)

    # ValueError: greater than limit
    limit = 4
    wrong_values = [5, 6]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must be smaller than'):
            _validate_int('vw', wv, limit=limit)

    # ValueError: wrong parity
    wrong_values = [1, 3]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must be even'):
            _validate_int('vw', wv, parity='even')
    wrong_values = [2, 4]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must be odd'):
            _validate_int('vw', wv, parity='odd')
    # Type: wrong parity string
    wrong_values = [2, 4]
    for wv in wrong_values:
        with pytest.raises(TypeError, match='unknown parity value'):
            _validate_int('vw', wv, parity='asdf')


def test_validate_float():

    # TypeError: not a float (NOTE: np.nan, np.inf are actually floats!)
    wrong_values = ['1', True, 1, 3j, [2.3], None, np.asarray([0.5])]
    for wv in wrong_values:
        with pytest.raises(TypeError, match='must be an instance of float'):
            _validate_float('vw', wv)

    # ValueError: not positive
    wrong_values = [-1.1, -3.14, -1e-100]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must be positive'):
            _validate_float('vw', wv, positive=True)

    # ValueError: greater than limit
    limit = 5.
    wrong_values = [5+1e-5, 6., np.inf]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must be smaller than'):
            _validate_float('vw', wv, limit=limit)


def test_validate_number():

    # TypeError: not a float, int or ndarray
    wrong_values = ['1', True, 3j, [2.3], None]
    for wv in wrong_values:
        with pytest.raises(TypeError, match='must be an integer, float or 1-D ndarray'):
            _validate_number('vw', wv)

    # ValueError: not norm
    wrong_values = [1.1, -0.5, 2, np.asarray([0.5, 1.1])]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must be in the interval'):
            _validate_number('vw', wv, norm=True)

    # ValueError: not positive
    wrong_values = [-1.1, -0.5, -2, np.asarray([-0.5, 1.1])]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must be positive'):
            _validate_number('vw', wv, positive=True)

    # ValueError: outside limits
    limit = [-1.5, 0.5]
    wrong_values = [-2, -1.6, 1, 0.51, np.asarray([-3, 0]), np.asarray([-0.5, 1.])]
    for wv in wrong_values:
        with pytest.raises(ValueError):
            _validate_number('vw', wv, limit=limit)


def test_validate_list():

    # TypeError: not a list
    wrong_values = ['1', True, 3j, np.asarray([2.3]), None, np.nan, np.inf]
    for wv in wrong_values:
        with pytest.raises(TypeError, match='must be an instance of list'):
            _validate_list('vw', wv)

    # Value error: size mismatch
    size = 3
    wrong_values = [ [], [1], [None, None], [1,2,3,4]]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must have size'):
            _validate_list('vw', wv, size=size)


def test_validate_ndarray():

    # TypeError: not a ndarray
    wrong_values = ['1', True, 3j, [2.3], None, np.nan, np.inf]
    for wv in wrong_values:
        with pytest.raises(TypeError, match='must be an instance of ndarray'):
            _validate_ndarray('vw', wv)

def test_validate_ndarray_1D():

    # TypeError: not a ndarray
    wrong_values = ['1', True, 3j, [2.3], None, np.nan, np.inf]
    for wv in wrong_values:
        with pytest.raises(TypeError, match='must be an instance of ndarray'):
            _validate_ndarray_1D('vw', wv)

    # ValueError: not 1D
    wrong_values = [np.empty(()), np.empty((1,2)), np.empty((1,2,3))]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must be 1D'):
            _validate_ndarray_1D('vw', wv)

    # ValueError: not given size
    size = 3
    wrong_values = [np.empty((1)), np.empty((2)), np.empty((4))]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must have size'):
            _validate_ndarray_1D('vw', wv, size=size)

    # ValueError: not norm
    wrong_values = [np.ones((3))*2, np.ones((3))*-1, np.asarray([0., 1., 2.])]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must be in the interval'):
            _validate_ndarray_1D('vw', wv, norm=True)

    # ValueError: not positive
    wrong_values = [np.ones((3))*-1, np.asarray([0., -1., 2.])]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must be positive'):
            _validate_ndarray_1D('vw', wv, positive=True)

    # ValueError: outside limits
    limit = [-1.5, 0.5]
    wrong_values = [np.ones((3)), np.asarray([0., -1., -2.])]
    for wv in wrong_values:
        with pytest.raises(ValueError):
            _validate_ndarray_1D('vw', wv, limit=limit)

    # ValueError: wrong dtype
    dtype = np.dtype(float)
    wrong_values = [np.ones((3), dtype='int'), np.asarray([0., -1., 2.], dtype='O')]
    for wv in wrong_values:
        with pytest.raises(TypeError, match='dtype must be'):
            _validate_ndarray_1D('vw', wv, dtype=dtype)


def test_validate_ndarray_2D():

    # TypeError: not a ndarray
    wrong_values = ['1', True, 3j, [2.3], None, np.nan, np.inf]
    for wv in wrong_values:
        with pytest.raises(TypeError, match='must be an instance of ndarray'):
            _validate_ndarray_2D('vw', wv)

    # ValueError: not 2D
    wrong_values = [np.empty(()), np.empty((1)), np.empty((1,2,3))]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must be 2D'):
            _validate_ndarray_2D('vw', wv)

    # ValueError: shape mismatch
    shape0 = 3
    wrong_values = [np.empty((1,3)), np.empty((2,3)), np.empty((4,3))]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must have dimension 0='):
            _validate_ndarray_2D('vw', wv, shape0=shape0)
    shape1 = 1
    wrong_values = [np.empty((1,3)), np.empty((2,3)), np.empty((4,3))]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must have dimension 1='):
            _validate_ndarray_2D('vw', wv, shape1=shape1)

    # ValueError: not norm
    wrong_values = [np.ones((3,2))*2, np.ones((3,2))*-1, np.asarray([[0., 1., 2.], [3., 3., 3.]])]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must be in the interval'):
            _validate_ndarray_2D('vw', wv, norm=True)

    # ValueError: not positive
    wrong_values = [np.ones((3,2))*-1, np.asarray([[0., -1., 2.], [3., 3., 3.]])]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must be positive'):
            _validate_ndarray_2D('vw', wv, positive=True)

    # ValueError: outside limits
    limit = [-1.5, 0.5]
    wrong_values = [np.ones((3,2)), np.asarray([[0., -1., 2.], [3., 3., -3.]])]
    for wv in wrong_values:
        with pytest.raises(ValueError):
            _validate_ndarray_2D('vw', wv, limit=limit)

    # ValueError: wrong dtype
    dtype = np.dtype(float)
    wrong_values = [np.ones((3,2), dtype='int'), np.asarray([[0., -1., 2.], [3., 3., 3.]], dtype='O')]
    for wv in wrong_values:
        with pytest.raises(TypeError, match='dtype must be'):
            _validate_ndarray_2D('vw', wv, dtype=dtype)


def test_validate_ndarray_3D():

    # TypeError: not a ndarray
    wrong_values = ['1', True, 3j, [2.3], None, np.nan, np.inf]
    for wv in wrong_values:
        with pytest.raises(TypeError, match='must be an instance of ndarray'):
            _validate_ndarray_3D('vw', wv)

    # ValueError: not 3D
    wrong_values = [np.empty(()), np.empty((1)), np.empty((1,2)), np.empty((1,2,3,4))]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must be 3D'):
            _validate_ndarray_3D('vw', wv)

    # ValueError: shape mismatch
    shape0 = 3
    wrong_values = [np.empty((1,2,3)), np.empty((2,2,2)), np.empty((4,3,1))]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must have dimension 0='):
            _validate_ndarray_3D('vw', wv, shape0=shape0)
    shape1 = 1
    wrong_values = [np.empty((1,2,3)), np.empty((2,2,2)), np.empty((4,3,1))]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must have dimension 1='):
            _validate_ndarray_3D('vw', wv, shape1=shape1)
    shape2 = 4
    wrong_values = [np.empty((1, 2, 3)), np.empty((2, 2, 2)), np.empty((4, 3, 1))]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must have dimension 2='):
            _validate_ndarray_3D('vw', wv, shape2=shape2)

    # ValueError: not norm
    wrong_values = [np.ones((3,2,1))*2, np.ones((3,2,1))*-1]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must be in the interval'):
            _validate_ndarray_3D('vw', wv, norm=True)

    # ValueError: not positive
    wrong_values = [np.ones((3,2,1))*-1]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must be positive'):
            _validate_ndarray_3D('vw', wv, positive=True)

    # ValueError: outside limits
    limit = [-1.5, 0.5]
    wrong_values = [np.ones((3,2,1)), np.ones((3,2,1))*-2]
    for wv in wrong_values:
        with pytest.raises(ValueError):
            _validate_ndarray_3D('vw', wv, limit=limit)

    # ValueError: wrong dtype
    dtype = np.dtype(float)
    wrong_values = [np.ones((3,2,1), dtype='int'), np.empty((1,2,3), dtype='O')]
    for wv in wrong_values:
        with pytest.raises(TypeError, match='dtype must be'):
            _validate_ndarray_3D('vw', wv, dtype=dtype)


def test_validate_ndarray_4D():

    # TypeError: not a ndarray
    wrong_values = ['1', True, 3j, [2.3], None, np.nan, np.inf]
    for wv in wrong_values:
        with pytest.raises(TypeError, match='must be an instance of ndarray'):
            _validate_ndarray_4D('vw', wv)

    # ValueError: not 4D
    wrong_values = [np.empty(()), np.empty((1)), np.empty((1,2)), np.empty((1,2,3)), np.empty((1,2,3,4,5))]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must be 4D'):
            _validate_ndarray_4D('vw', wv)

    # ValueError: shape mismatch
    shape0 = 3
    wrong_values = [np.empty((1,2,3,4)), np.empty((2,2,2,2)), np.empty((4,3,1,1))]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must have dimension 0='):
            _validate_ndarray_4D('vw', wv, shape0=shape0)
    shape1 = 1
    wrong_values = [np.empty((1,2,3,4)), np.empty((2,2,2,2)), np.empty((4,3,1,1))]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must have dimension 1='):
            _validate_ndarray_4D('vw', wv, shape1=shape1)
    shape2 = 4
    wrong_values = [np.empty((1, 2, 3,4)), np.empty((2, 2, 2, 2)), np.empty((4, 3, 1, 1))]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must have dimension 2='):
            _validate_ndarray_4D('vw', wv, shape2=shape2)
    shape3 = 3
    wrong_values = [np.empty((1, 2, 3, 4)), np.empty((2, 2, 2, 2)), np.empty((4, 3, 1, 1))]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must have dimension 3='):
            _validate_ndarray_4D('vw', wv, shape3=shape3)

    # ValueError: not norm
    wrong_values = [np.ones((3,2,1,1))*2, np.ones((3,2,1,1))*-1]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must be in the interval'):
            _validate_ndarray_4D('vw', wv, norm=True)

    # ValueError: not positive
    wrong_values = [np.ones((3,2,1,1))*-1]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must be positive'):
            _validate_ndarray_4D('vw', wv, positive=True)

    # ValueError: outside limits
    limit = [-1.5, 0.5]
    wrong_values = [np.ones((3,2,1,1)), np.ones((3,2,1,1))*-2]
    for wv in wrong_values:
        with pytest.raises(ValueError):
            _validate_ndarray_4D('vw', wv, limit=limit)

    # ValueError: wrong dtype
    dtype = np.dtype(float)
    wrong_values = [np.ones((3,2,1,1), dtype='int'), np.zeros((1,2,3,4), dtype='O')]
    for wv in wrong_values:
        with pytest.raises(TypeError, match='dtype must be'):
            _validate_ndarray_4D('vw', wv, dtype=dtype)


def test_validate_string():

    # TypeError: not a list
    wrong_values = [1, 3.14, True, 3j, np.asarray([2.3]), None, np.nan, np.inf]
    for wv in wrong_values:
        with pytest.raises(TypeError, match='must be an instance of str'):
            _validate_string('vw', wv)

    # Value error: not a given choice
    choices = ['foo', 'bar', 'ambisonic']
    wrong_values = ['FOO', 'ambisonics', 'random_string']
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must be one of the following'):
            _validate_string('vw', wv, choices=choices)


def test_validate_echogram():

    # TypeError: not an Echogram
    wrong_values = ['1', 1, 3.14, True, 3j, np.asarray([2.3]), None, np.nan, np.inf, QuantisedEchogram(None, None, None)]
    for wv in wrong_values:
        with pytest.raises(TypeError, match='echogram must be an instance of Echogram'):
            _validate_echogram(wv)

    # ValueError: shape mismatch
    with pytest.raises(ValueError, match='echogram shape mismatch'):
        l = 100
        value = np.ones((l, 1))
        time = np.ones((l))
        order = np.ones((l, C), dtype='int')
        coords = np.ones((l+1, C))
        wv = Echogram(value, time, order, coords)
        _validate_echogram(wv)


def test_validate_quantised_echogram():

    # TypeError: not an Echogram
    wrong_values = ['1', 1, 3.14, True, 3j, np.asarray([2.3]), None, np.nan, np.inf, Echogram(None, None, None, None)]
    for wv in wrong_values:
        with pytest.raises(TypeError, match='quantised echogram must be an instance of QuantisedEchogram'):
            _validate_quantised_echogram(wv)

    # ValueError: shape mismatch
    with pytest.raises(ValueError, match='quantised echogram shape mismatch'):
        l = 100
        value = np.ones((l, 1))
        time = np.ones((l+1))
        isActive = True
        wv = QuantisedEchogram(value, time, isActive)
        _validate_quantised_echogram(wv)


def test_validate_echogram_array():

    # TypeError: not a ndarray
    wrong_values = ['1', 1, 3.14, True, 3j, None, np.nan, np.inf, Echogram(None, None, None, None), QuantisedEchogram(None, None, None)]
    for wv in wrong_values:
        with pytest.raises(TypeError, match='Echogram array must be an instance of ndarray'):
            _validate_echogram_array(wv)

    # ValueError: dim not 2 or 3
    wrong_values = [np.ones(()), np.ones((1)), np.ones((1,2,3,4))]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='Echogram array must be 2D or 3D'):
            _validate_echogram_array(wv)

    # ValueError: shape mismatch
    shape0 = 3
    wrong_values = [np.empty((1,2,3)), np.empty((2,2,2)), np.empty((4,3,1))]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='Echogram array must have dimension 0='):
            _validate_echogram_array(wv, shape0=shape0)
    shape1 = 1
    wrong_values = [np.empty((1,2,3)), np.empty((2,2,2)), np.empty((4,3,1))]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='Echogram array must have dimension 1='):
            _validate_echogram_array(wv, shape1=shape1)
    shape2 = 4
    wrong_values = [np.empty((1, 2, 3)), np.empty((2, 2, 2,)), np.empty((4, 3, 1))]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='Echogram array must have dimension 2='):
            _validate_echogram_array(wv, shape2=shape2)

    # TypeError: wrong dtype
    dtype = np.dtype(float)
    wrong_values = [np.zeros((1,2,3), dtype='float'), np.ones((2,2,2), dtype='int'), np.empty((4,3,1), dtype='O')]
    for wv in wrong_values:
        with pytest.raises(TypeError, match='Echogram array dtype must be Echogram'):
            _validate_echogram_array(wv)

    # ValueError: invalid echograms
    with pytest.raises(ValueError, match='echogram shape mismatch'):
        l = 100
        value = np.ones((l, 1))
        time = np.ones((l))
        order = np.ones((l, C), dtype='int')
        coords = np.ones((l+1, C))
        e = Echogram(value, time, order, coords)
        wv = np.asarray([[e, e], [e, e]])
        _validate_echogram_array(wv)

def test_validate_quantised_echogram_array():

    # TypeError: not a ndarray
    wrong_values = ['1', 1, 3.14, True, 3j, None, np.nan, np.inf, Echogram(None, None, None, None), QuantisedEchogram(None, None, None)]
    for wv in wrong_values:
        with pytest.raises(TypeError, match='Quantised Echogram array must be an instance of ndarray'):
            _validate_quantised_echogram_array(wv)

    # ValueError: dim not 1
    wrong_values = [np.ones(()), np.ones((1,2)), np.ones((1,2,3))]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='Quantised Echogram array must be 1D'):
            _validate_quantised_echogram_array(wv)

    # ValueError: size mismatch
    size = 3
    wrong_values = [np.empty((1)), np.empty((2)), np.empty((4))]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='Quantised Echogram array must have size='):
            _validate_quantised_echogram_array(wv, size=size)

    # TypeError: wrong dtype
    dtype = np.dtype(float)
    wrong_values = [np.zeros((1), dtype='float'), np.ones((2), dtype='int'), np.empty((4), dtype='O')]
    for wv in wrong_values:
        with pytest.raises(TypeError, match='Quantised Echogram array dtype must be QuantisedEchogram'):
            _validate_quantised_echogram_array(wv)

    # ValueError: invalid quantised echograms
    with pytest.raises(ValueError, match='quantised echogram shape mismatch'):
        l = 100
        value = np.ones((l, 1))
        time = np.ones((l+1))
        isActive = True
        e = QuantisedEchogram(value, time, isActive)
        wv = np.asarray([e, e])
        _validate_quantised_echogram_array(wv)
