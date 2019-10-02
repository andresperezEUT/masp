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
#   @file   test_sph_array_characteristics.py
#   @author Andrés Pérez-López
#   @date   31/07/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


import random
from masp.tests.convenience_test_methods import *

def test_sph_array_noise():
    num_tests = 10
    params = {
        'R':
        [np.random.random() for i in range(num_tests)],
        'Nmic':
        [np.random.randint(1,10) for i in range(num_tests)],
        'maxN':
        [np.random.randint(1,10) for i in range(num_tests)],
        'arrayType':
        [random.choice(['open', 'rigid']) for i in range(num_tests)],
        'f':
        [(np.arange(10) * np.random.randint(100000) + 100).tolist() for i in range(num_tests)],
    }
    for t in range(num_tests):
        p = get_parameters(params, t)
        numeric_assert("sphArrayNoise",
                       "sph_array_noise",
                       *p,
                       nargout=2,
                       namespace='sap')


def test_sph_array_noise_threshold():
    num_tests = 10
    params = {
        'R':
        [np.random.random() for i in range(num_tests)],
        'Nmic':
        [np.random.randint(1,10) for i in range(num_tests)],
        'maxG_db':
        [np.random.random()*20-10 for i in range(num_tests)],
        'maxN':
        [np.random.randint(1,10) for i in range(num_tests)],
        'arrayType':
        [random.choice(['open', 'rigid', 'directional']) for i in range(num_tests)],
        'dirCoef':
        [np.random.random() for i in range(num_tests)]
    }
    for t in range(num_tests):
        p = get_parameters(params, t)
        numeric_assert("sphArrayNoiseThreshold",
                       "sph_array_noise_threshold",
                       *p,
                       nargout=1,
                       namespace='sap')


def test_sph_array_alias_lim():
    num_tests = 10
    nMics = [np.random.randint(1, 10) for i in range(num_tests)]
    params = {
        'R':
        [np.random.random() for i in range(num_tests)],
        'Nmic':
        nMics,
        'maxN':
        [np.random.randint(1,10) for i in range(num_tests)],
        'mic_dirs_rad':
        [(np.random.rand(np.random.randint(1,10),C-1)*[2*np.pi, np.pi]).tolist() for i in range(num_tests)],
        'mic_weights':
        [None for i in range(num_tests)]  # TODO not implemented
    }
    for t in range(num_tests):
        p = get_parameters(params, t)
        numeric_assert("sphArrayAliasLim",
                       "sph_array_alias_lim",
                       *p,
                       nargout=1,
                       namespace='sap')
