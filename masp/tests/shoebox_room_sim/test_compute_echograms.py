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
#   @file   test_compute_echograms.py
#   @author Andrés Pérez-López
#   @date   14/08/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from masp.tests.convenience_test_methods import *
from masp.utils import C


def test_compute_echograms_mic():
    num_tests = 5
    nBands = [np.random.randint(1, 5) for i in range(num_tests)]
    nSrc = [np.random.randint(1, 5) for i in range(num_tests)]
    nRec = [np.random.randint(1, 5) for i in range(num_tests)]
    params = [
        # room
        [(np.random.random(C) * 5 + 5).tolist() for i in range(num_tests)],
        # source
        [(np.random.random((nSrc[i],C)) * 5).tolist() for i in range(num_tests)],
        # receiver
        [(np.random.random((nRec[i],C)) * 5).tolist() for i in range(num_tests)],
        # abs_wall
        [np.random.random((nBands[i], 2*C)).tolist() for i in range(num_tests)],
        # limits
        [(np.random.random(nBands[i]) + 0.1).tolist() for i in range(num_tests)],
        # mic specs
        [generate_random_mic_specs(nRec[i]) for i in range(num_tests)],
    ]
    num_params = len(params)
    for t in range(num_tests):
        print('-----------------------------------------------')
        print('                  t='+str(t))
        p = []
        for p_idx in range(num_params):
            p.append(params[p_idx][t])
            print(params[p_idx][t])
        print('-----------------------------------------------')
        numeric_assert("compute_echograms_mics",
                       "compute_echograms_mic",
                       *p,
                       write_file=True,
                       namespace='srs')
