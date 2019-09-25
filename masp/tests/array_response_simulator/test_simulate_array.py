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
#   @file   test_sph_functions.py
#   @author Andrés Pérez-López
#   @date   22/08/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import pytest
from masp.tests.convenience_test_methods import *
from masp.utils import C
import random

def test_simulate_sph_array():
    num_tests = 10
    nMics = [np.random.randint(1, 10) for i in range(num_tests)]
    nDOAs = [np.random.randint(1, 10) for i in range(num_tests)]
    params = {
        'N_filt': # even
        [np.random.randint(10, 100)*2 for i in range(num_tests)],
        'mic_dirs_rad':
        [(np.random.rand(nMics[i], C-1)*[2*np.pi, np.pi]).tolist() for i in range(num_tests)],
        'src_dirs_rad':
        [(np.random.rand(nDOAs[i], C-1)*[2*np.pi, np.pi]).tolist() for i in range(num_tests)],
        'arrayType':
        [random.choice(['open', 'rigid', 'directional']) for i in range(num_tests)],
        'R':
        [np.random.random() for i in range(num_tests)],
        'N_order':
        [np.random.randint(1, 10) for i in range(num_tests)],
        'fs':
        [np.random.randint(100000) + 100 for i in range(num_tests)],
        'dirCoef':
        [np.random.random() for i in range(num_tests)],
    }
    for t in range(num_tests):
        p = get_parameters(params, t)
        numeric_assert("simulateSphArray",
                       "simulate_sph_array",
                       *p,
                       nargout=2,
                       namespace='ars')

    # ValueError: directional array, but no dirCoef defined
    with pytest.raises(ValueError, match='dirCoef must be defined in the directional case'):
        N_filt = p[0]
        mic_dirs_rad = np.asarray(p[1])
        src_dirs_rad = np.asarray(p[2])
        arrayType = 'directional'
        R = p[4]
        N_order = p[5]
        fs = p[6]
        dirCoef = None
        masp.ars.simulate_sph_array(N_filt, mic_dirs_rad, src_dirs_rad, arrayType, R, N_order, fs, dirCoef)


def test_simulate_cyl_array():
    num_tests = 10
    nMics = [np.random.randint(1, 10) for i in range(num_tests)]
    nDOAs = [np.random.randint(1, 10) for i in range(num_tests)]
    params = {
        'N_filt': # even
        [np.random.randint(10, 100)*2 for i in range(num_tests)],
        'mic_dirs_rad':
        [(np.random.rand(nMics[i])*[2*np.pi]).tolist() for i in range(num_tests)],
        'src_dirs_rad':
        [(np.random.rand(nDOAs[i])*[2*np.pi]).tolist() for i in range(num_tests)],
        'arrayType':
        [random.choice(['open', 'rigid']) for i in range(num_tests)],
        'R':
        [np.random.random() for i in range(num_tests)],
        'N_order':
        [np.random.randint(1, 10) for i in range(num_tests)],
        'fs':
        [np.random.randint(100000) + 100 for i in range(num_tests)],
    }
    for t in range(num_tests):
        p = get_parameters(params, t)
        numeric_assert("simulateCylArray",
                       "simulate_cyl_array",
                       *p,
                       nargout=2,
                       namespace='ars')


def test_get_array_response():
    num_tests = 10
    nMics = [np.random.randint(1, 10) for i in range(num_tests)]
    nDOAs = [np.random.randint(1, 10) for i in range(num_tests)]

    params = {
        'src_dirs':
        [(np.random.rand(nDOAs[i],C)*[1, 1, 0]).tolist() for i in range(num_tests)],
        'mic_pos':
        [(np.random.rand(nMics[i],C) * [1, 1, 0]).tolist() for i in range(num_tests)],
        'N_filt':
        [np.random.randint(10, 100) * 2 for i in range(num_tests)],
        'fs':
        [np.random.randint(100000) + 100 for i in range(num_tests)],
        'mic_dirs':
        [random.choice([
            None,
            (np.random.rand(C)).tolist(),
            (np.random.rand(nMics[i], C)).tolist(),
        ]) for i in range(num_tests)],
        'fDir_handle':
        # TODO: how to pass lambda expressions to Matlab...?
        [None for i in range(num_tests)],

    }
    for t in range(num_tests):
        p = get_parameters(params, t)
        numeric_assert("getArrayResponse",
                       "get_array_response",
                       *p,
                       nargout=2,
                       namespace='ars')
