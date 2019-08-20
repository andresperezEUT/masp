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

from masp.tests.convenience_test_methods import *
import random


def test_render_rirs():
    num_tests = 10
    params = {
        'echogram':
        [generate_random_echogram() for i in range(num_tests)],
        'endtime':
        [np.random.rand()+0.01 for i in range(num_tests)],
        # [0.01 for i in range(num_tests)],
        'fs':
        [np.random.randint(100000)+100 for i in range(num_tests)],
        # [2000 for i in range(num_tests)],
        'fractional': 
        [random.choice([True, False]) for i in range(num_tests)]
        # [True for i in range(num_tests)]
    }
    for t in range(num_tests):
        p = get_parameters(params, t)
        numeric_assert("render_rir",
                       "render_rirs",
                       *p,
                       nargout=1,
                       namespace='srs')


# test_render_rirs()


#
# def test_render_rirs_mic():
#
#     raise NotImplementedError
#
#     num_tests = 5
#     nSrc = [np.random.randint(1, 10) for i in range(num_tests)]
#     nRec = [np.random.randint(1, 10) for i in range(num_tests)]
#     max_sh_order = 10
#     params = {
#         'echograms':
#             [generate_random_echogram_array(nSrc[i], nRec[i]) for i in range(num_tests)],
#         'sh_order':
#             [random.choice([np.random.randint(max_sh_order),
#                             (np.random.rand(nRec[i]) * max_sh_order).round().astype(int).tolist()]) for i in
#              range(num_tests)],
#     }
#     for t in range(num_tests):
#         p = get_parameters(params, t)
#         numeric_assert("render_mic_rirs_test",
#                        "render_rirs_mic",
#                        *p,
#                        nargout=1,
#                        namespace='srs')
