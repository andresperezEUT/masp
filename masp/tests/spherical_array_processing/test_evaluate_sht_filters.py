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
#   @file   test_evaluate_sht_filters.py
#   @author Andrés Pérez-López
#   @date   31/07/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import numpy as np
import random
from masp.tests.convenience_test_methods import *


def test_evaluate_sht_filters():
    num_tests = 10
    sh_order = [np.random.randint(10) for i in range(num_tests)]
    n_sh = [int(np.power(sh_order[i]+1, 2)) for i in range(num_tests)]
    nMics = [np.random.randint(2,10) for i in range(num_tests)]
    nBins = [np.random.randint(1,10) for i in range(num_tests)]
    nGrid = [np.random.randint(2,10) for i in range(num_tests)]
    print('sh_order', sh_order)
    print('n_sh', n_sh)
    print('nMics', nMics)
    print('nBins', nBins)
    print('nGrid', nGrid)
    params = {
        'M_mic2sh':
        [np.random.rand(n_sh[i], nMics[i], nBins[i]).tolist() for i in range(num_tests)],
        'H_array':  # complex
        [((np.random.rand(nBins[i],nMics[i],nGrid[i])+np.random.rand(nBins[i],nMics[i],nGrid[i])*1j)*2-1).tolist() for i in range(num_tests)],
        'fs':
        [np.random.randint(100000) + 100 for i in range(num_tests)],
        'Y_grid':
        [np.random.rand(nGrid[i], n_sh[i]).tolist() for i in range(num_tests)],
        'w_grid':
        [random.choice([None, np.random.rand(nGrid[i]).tolist()]) for i in range(num_tests)],
    }
    for t in range(num_tests):
        p = get_parameters(params, t)
        numeric_assert("evaluateSHTfilters",
                       "evaluate_sht_filters",
                       *p,
                       nargout=3,
                       namespace='sap')
