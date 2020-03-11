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
#   @file   test_render_rirs.py
#   @author Andrés Pérez-López
#   @date   30/07/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
from asma.shoebox_room_sim import quantise_echogram
from asma.tests.convenience_test_methods import *
import random

def test_render_rirs_array():
    num_tests = 5
    nSrc = [np.random.randint(1, 5) for i in range(num_tests)]
    nRec = [np.random.randint(1, 5) for i in range(num_tests)]
    nBands = [np.random.randint(1, 5) for i in range(num_tests)]
    band_centerfreqs = [[np.random.randint(30, 100)] for i in range(num_tests)]
    [[band_centerfreqs[i].append(2 * band_centerfreqs[i][b]) for b in range(nBands[i] - 1)] for i in range(num_tests)]
    fs = [int(band_centerfreqs[i][-1]*2.1) for i in range(num_tests)] # avoid fs/2 smaller than last centerband

    num_grid_points = [[np.random.randint(1,10) for r in range(nRec[i])] for i in range(num_tests)]
    num_mics = [[np.random.randint(1,10) for r in range(nRec[i])] for i in range(num_tests)]

    # Explicit conversion of list to tuples: 1D matlab cell arrays
    grids = []
    array_irs = []
    for i in range(num_tests):
        grid = [(np.random.rand(num_grid_points[i][r], C-1) * [2*np.pi, np.pi]).tolist() for r in range(nRec[i])]
        grids.append(tuple(grid))
        array_ir = [(np.random.rand(256, num_mics[i][r], num_grid_points[i][r])).tolist() for r in range(nRec[i])]
        array_irs.append(tuple(array_ir))

    params = {
        'echograms':
        [generate_random_echogram_array(nSrc[i], nRec[i], nBands[i]) for i in range(num_tests)],
        'band_centerfreqs':
        band_centerfreqs,
        'fs':
        fs,
        'grids':
        grids,
        'array_irs':
        array_irs,
    }
    for t in range(num_tests):
        p = get_parameters(params, t)
        numeric_assert("render_array_rirs",
                       "render_rirs_array",
                       *p,
                       nargout=1,
                       namespace='srs')


def test_render_rirs_mic():
    num_tests = 5
    nSrc = [np.random.randint(1, 5) for i in range(num_tests)]
    nRec = [np.random.randint(1, 5) for i in range(num_tests)]
    nBands = [np.random.randint(1, 5) for i in range(num_tests)]
    band_centerfreqs = [[np.random.randint(30, 100)] for i in range(num_tests)]
    [[band_centerfreqs[i].append(2 * band_centerfreqs[i][b]) for b in range(nBands[i] - 1)] for i in range(num_tests)]
    fs = [int(band_centerfreqs[i][-1]*2.1) for i in range(num_tests)] # avoid fs/2 smaller than last centerband

    params = {
        'echograms':
        [generate_random_echogram_array(nSrc[i], nRec[i], nBands[i]) for i in range(num_tests)],
        'band_centerfreqs':
        band_centerfreqs,
        'fs':
        fs
    }
    for t in range(num_tests):
        p = get_parameters(params, t)
        numeric_assert("render_mic_rirs",
                       "render_rirs_mic",
                       *p,
                       nargout=1,
                       namespace='srs')


def test_render_rirs_sh():
    num_tests = 5
    nSrc = [np.random.randint(1, 5) for i in range(num_tests)]
    nRec = [np.random.randint(1, 5) for i in range(num_tests)]
    nBands = [np.random.randint(1, 5) for i in range(num_tests)]
    band_centerfreqs = [[np.random.randint(30, 100)] for i in range(num_tests)]
    [[band_centerfreqs[i].append(2 * band_centerfreqs[i][b]) for b in range(nBands[i] - 1)] for i in range(num_tests)]
    fs = [int(band_centerfreqs[i][-1]*2.1) for i in range(num_tests)] # avoid fs/2 smaller than last centerband

    params = {
        'echograms':
        [generate_random_echogram_array_sh(nSrc[i], nRec[i], nBands[i]) for i in range(num_tests)],
        'band_centerfreqs':
        band_centerfreqs,
        'fs':
        fs
    }
    for t in range(num_tests):
        p = get_parameters(params, t)
        numeric_assert("render_sh_rirs",
                       "render_rirs_sh",
                       *p,
                       nargout=1,
                       namespace='srs')


def test_render_rirs():
    num_tests = 5
    params = {
        'echogram':
        [generate_random_echogram() for i in range(num_tests)],
        'endtime':
        [np.random.rand()+0.01 for i in range(num_tests)],
        'fs':
        [np.random.randint(100000)+100 for i in range(num_tests)],
        'fractional': 
        [random.choice([True, False]) for i in range(num_tests)]
    }
    for t in range(num_tests):
        p = get_parameters(params, t)
        numeric_assert("render_rir",
                       "render_rirs",
                       *p,
                       nargout=1,
                       namespace='srs')

def test_render_quantised():
    num_tests = 10
    echograms = [generate_random_echogram() for i in range(num_tests)]
    nGrids = [np.random.randint(1, 10) for i in range(num_tests)]
    echo2GridMaps = [np.random.randint(0, 100, echograms[i].time.size + 1) for i in range(num_tests)]

    params = {
        'qechogram':
        [quantise_echogram(echograms[i], nGrids[i], echo2GridMaps[i]) for i in range(num_tests)],
        'endtime':
        [np.random.rand() + 0.01 for i in range(num_tests)],
        'fs':
        [np.random.randint(100000) + 100 for i in range(num_tests)],
        'fractional':
        [random.choice([True, False]) for i in range(num_tests)]
    }
    for t in range(num_tests):
        p = get_parameters(params, t)
        numeric_assert("render_quantized",
                       "render_quantised",
                       *p,
                       nargout=2,
                       namespace='srs')


def test_filter_rirs():
    num_tests = 5
    nBands = [np.random.randint(1, 10) for i in range(num_tests)]
    band_centerfreqs = [[np.random.randint(30,100)] for i in range(num_tests)]
    [[band_centerfreqs[i].append(2 * band_centerfreqs[i][b]) for b in range(nBands[i] - 1)] for i in range(num_tests)]
    nFrames = [np.random.randint(1, 10)*100 for i in range(num_tests)]
    fs = [int(band_centerfreqs[i][-1]*2.1) for i in range(num_tests)] # avoid fs/2 smaller than last centerband
    params = {
        'rir':
        [(np.random.rand(nFrames[i], nBands[i])*2-1).tolist() for i in range(num_tests)],
        'f_center':
        band_centerfreqs,
        'fs':
        fs,
    }
    for t in range(num_tests):
        p = get_parameters(params, t)
        numeric_assert("filter_rir",
                       "filter_rirs",
                       *p,
                       nargout=1,
                       namespace='srs')
