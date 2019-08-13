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
#   @file   test_acoustics.py
#   @author Andrés Pérez-López
#   @date   31/07/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


from masp.tests.convenience_test_methods import *
import random


def test_find_abs_coeffs_from_rt():
    num_tests = 10
    params = [
        # room
        [(np.random.random(3) * 10).tolist() for i in range(num_tests)],
        # rt60
        [(np.random.random(6)).tolist() for i in range(num_tests)],
        # abs_wall_ratios
        [random.choice([None, (np.random.random(6) * 10).tolist()]) for i in range(num_tests)],
        # [None for i in range(num_tests)],
    ]
    num_params = len(params)
    for t in range(num_tests):
        p = []
        for p_idx in range(num_params):
            p.append(params[p_idx][t])
        numeric_assert("findAbsCoeffsFromRT",
                       "find_abs_coeffs_from_rt",
                       *p,
                       nargout=2,
                       namespace='srs')


def test_get_rt_sabine():
    num_tests = 10
    params = [
        # alpha
        [np.random.rand() for i in range(num_tests)],
        # room
        [(np.random.random(3) * 10).tolist() for i in range(num_tests)],
        # abs_wall_ratios
        [(np.random.random(6) * 10).tolist() for i in range(num_tests)],
    ]
    num_params = len(params)
    for t in range(num_tests):
        p = []
        for p_idx in range(num_params):
            p.append(params[p_idx][t])
        numeric_assert("getRTsabine",
                       "get_rt_sabine",
                       *p,
                       nargout=1,
                       namespace='srs')


def test_room_stats():
    num_tests = 10
    params = [
        # room
        [(np.random.random(3) * 10).tolist() for i in range(num_tests)],
        # abs_wall
        [np.random.random((np.random.randint(1,10),6)).tolist() for i in range(num_tests)],
    ]
    num_params = len(params)
    for t in range(num_tests):
        p = []
        for p_idx in range(num_params):
            p.append(params[p_idx][t])
        numeric_assert("room_stats",
                       "room_stats",
                       *p,
                       nargout=3,
                       namespace='srs')

