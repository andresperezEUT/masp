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
#   @file   test_array_sht_filters.py
#   @author Andrés Pérez-López
#   @date   03/10/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import random
from masp.tests.convenience_test_methods import *


def test_array_sht_filters_theory_radInverse():
    num_tests = 10
    params = {
        'R':
        [np.random.random() for i in range(num_tests)],
        'Nmic':
        [np.random.randint(1,10) for i in range(num_tests)],
        'order_sht':
        [np.random.randint(10) for i in range(num_tests)],
        'Lfilt':  # even
        [np.random.randint(10, 100) * 2 for i in range(num_tests)],
        'fs':
        [np.random.randint(100000) + 100 for i in range(num_tests)],
        'amp_threshold':
        [(np.random.random()*2-1)*20 for i in range(num_tests)],
    }
    for t in range(num_tests):
        p = get_parameters(params, t)
        numeric_assert("arraySHTfiltersTheory_radInverse",
                       "array_sht_filters_theory_radInverse",
                       *p,
                       nargout=2,
                       namespace='sap')

def test_array_sht_filters_theory_softLim():
    num_tests = 10
    params = {
        'R':
        [np.random.random() for i in range(num_tests)],
        'Nmic':
        [np.random.randint(1,10) for i in range(num_tests)],
        'order_sht':
        [np.random.randint(10) for i in range(num_tests)],
        'Lfilt':  # even
        [np.random.randint(10, 100) * 2 for i in range(num_tests)],
        'fs':
        [np.random.randint(100000) + 100 for i in range(num_tests)],
        'amp_threshold':
        [(np.random.random()*2-1)*20 for i in range(num_tests)],
    }
    for t in range(num_tests):
        p = get_parameters(params, t)
        numeric_assert("arraySHTfiltersTheory_softLim",
                       "array_sht_filters_theory_softLim",
                       *p,
                       nargout=2,
                       namespace='sap')

def test_array_sht_filters_theory_regLS():
    num_tests = 10
    nMics = [np.random.randint(1,10) for i in range(num_tests)]
    params = {
        'R':  # array_order should be at least 2, which means that for a fs_min=100, R >= 1.091802909610402
        [np.random.random()+1.1 for i in range(num_tests)],
        'mic_dirs_rad':
        [(np.random.rand(nMics[i], C - 1) * [2 * np.pi, np.pi]).tolist() for i in range(num_tests)],
        'order_sht':
        [np.random.randint(10) for i in range(num_tests)],
        'Lfilt':  # even
        [np.random.randint(10,100) * 2 for i in range(num_tests)],
        'fs':
        [np.random.randint(1000) + 100 for i in range(num_tests)],
        'amp_threshold':
        [(np.random.random()*2-1)*20 for i in range(num_tests)],
    }
    for t in range(num_tests):
        p = get_parameters(params, t)
        numeric_assert("arraySHTfiltersTheory_regLS",
                       "array_sht_filters_theory_regLS",
                       *p,
                       nargout=2,
                       namespace='sap',
                       atol=1e-2)  # For some reason related with ifft, TODO check

def test_array_sht_filters_measure_regLS():
    num_tests = 10
    nMics = [np.random.randint(2,10) for i in range(num_tests)]
    nGrid = [np.random.randint(2,10) for i in range(num_tests)]
    nFFT = [np.random.randint(10, 100) * 2 for i in range(num_tests)]
    nBins = [nFFT[i]//2+1 for i in range(num_tests)]
    print('nMics', nMics)
    print('nBins', nBins)
    print('nFFT', nFFT)
    print('nGrid', nGrid)
    params = {
        'H_array':  # complex
        [((np.random.rand(nBins[i],nMics[i],nGrid[i])+np.random.rand(nBins[i],nMics[i],nGrid[i])*1j)*2-1).tolist() for i in range(num_tests)],
        'order_sht':
        [np.random.randint(10) for i in range(num_tests)],
        'grid_dirs_rad':
        [(np.random.random((nGrid[i], C - 1)) * [2 * np.pi, np.pi]).tolist() for i in range(num_tests)],
        'w_grid':
        [np.random.rand(nGrid[i]).tolist() for i in range(num_tests)],
        'nFFT':  # even
        nFFT,
        'amp_threshold':
        [(np.random.random()*2-1)*20 for i in range(num_tests)],
    }
    for t in range(num_tests):
        p = get_parameters(params, t)
        numeric_assert("arraySHTfiltersMeas_regLS",
                       "array_sht_filters_measure_regLS",
                       *p,
                       nargout=2,
                       namespace='sap')

def test_array_sht_filters_measure_regLSHD():
    num_tests = 100
    nFFT = [np.random.randint(10, 100) * 2 for i in range(num_tests)]
    nBins = [nFFT[i]//2+1 for i in range(num_tests)]
    order_sht = [np.random.randint(1,3) for i in range(num_tests)]
    min_nMics = [np.power((order_sht[i]+1), 2) for i in range(num_tests)]
    min_nGrid = [np.power(2*(order_sht[i]+1), 2)for i in range(num_tests)]
    nMics = [np.random.randint(min_nMics[i], min_nMics[i]*2) for i in range(num_tests)]
    nGrid = [np.random.randint(min_nGrid[i], min_nGrid[i]*2) for i in range(num_tests)]
    print('nMics', nMics)
    print('nBins', nBins)
    print('nFFT ', nFFT)
    print('order_sht', order_sht)
    print('nGrid', nGrid)
    params = {
        'H_array':  # complex
        [((np.random.rand(nBins[i],nMics[i],nGrid[i])+np.random.rand(nBins[i],nMics[i],nGrid[i])*1j)*2-1).tolist() for i in range(num_tests)],
        'order_sht':
        order_sht,
        'grid_dirs_rad':
        [(np.random.random((nGrid[i], C - 1)) * [2 * np.pi, np.pi]).tolist() for i in range(num_tests)],
        'w_grid':
        [np.random.rand(nGrid[i]).tolist() for i in range(num_tests)],
        'nFFT':  # even
        nFFT,
        'amp_threshold':
        [(np.random.random()*2-1)*20 for i in range(num_tests)],
    }
    for t in range(num_tests):
        p = get_parameters(params, t, verbose=True)
        numeric_assert("arraySHTfiltersMeas_regLSHD",
                       "array_sht_filters_measure_regLSHD",
                       *p,
                       nargout=2,
                       namespace='sap')
