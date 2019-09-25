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
from masp.validate_data_types import _validate_boolean
from masp.validate_data_types import _validate_int
from masp.validate_data_types import _validate_float


def test_validate_boolean():

    # TypeError: not a boolean
    wrong_values = ['True', 1, 2.3, 1e4, 3j, [True], None, np.nan, np.inf]
    for wv in wrong_values:
        with pytest.raises(TypeError, match='must be an instance of bool'):
            _validate_boolean('vw', wv)

def test_validate_int():

    # TypeError: not an integer
    wrong_values = ['1', True, 2.3, 1e4, 3j, [1], None, np.nan, np.inf]
    for wv in wrong_values:
        with pytest.raises(TypeError, match='must be an instance of int'):
            _validate_int('vw', wv)

    # ValueError: not positive
    wrong_values = [-1, -3]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must be positive'):
            _validate_int('vw', wv, positive=True)

    # ValueError: greater than limit
    wrong_values = [5, 6]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must be smaller than'):
            _validate_int('vw', wv, limit=4)

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
    wrong_values = ['1', True, 1, 3j, [2.3], None]
    for wv in wrong_values:
        with pytest.raises(TypeError, match='must be an instance of float'):
            _validate_float('vw', wv)

    # ValueError: not positive
    wrong_values = [-1.1, -3.14, -1e-100]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must be positive'):
            _validate_float('vw', wv, positive=True)

    # ValueError: greater than limit
    wrong_values = [5+1e-5, 6., np.inf]
    for wv in wrong_values:
        with pytest.raises(ValueError, match='must be smaller than'):
            _validate_float('vw', wv, limit=5.)
