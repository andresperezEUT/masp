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
#   @file   test_quantise.py
#   @author Andrés Pérez-López
#   @date   09/09/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


import numpy as np
from masp.tests.convenience_test_methods import *
from masp.utils import C

def test_get_echo2gridMap():
    num_tests = 10
    params = {
        'echogram':
        [generate_random_echogram() for i in range(num_tests)],
        'grid_dirs_rad':
        [(np.random.random((np.random.randint(10, 100), C-1))*[2*np.pi, np.pi]).tolist() for i in range(num_tests)],
    }
    for t in range(num_tests):
        p = get_parameters(params, t)
        numeric_assert("get_echo2gridMap_test",
                       "get_echo2gridMap",
                       *p,
                       nargout=1,
                       namespace='srs')

def test_quantise_echogram():
    num_tests = 10
    echograms = [generate_random_echogram() for i in range(num_tests)]
    params = {
        'echogram':
        echograms,
        'nGrid':
        [np.random.randint(1,10) for i in range(num_tests)],
        'echo2gridMap':
        [(np.random.randint(0,100, echograms[i].time.size+1)).tolist() for i in range(num_tests)],
    }
    for t in range(num_tests):
        p = get_parameters(params, t)
        numeric_assert("quantise_echogram",
                       "quantise_echogram",
                       *p,
                       write_file=True,
                       namespace='srs')
