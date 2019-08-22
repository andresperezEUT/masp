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

from masp.tests.convenience_test_methods import *


def test_sph_besselj():
    num_tests = 10
    params = {
        'n':
        [np.random.randint(10) for i in range(num_tests)],
        'x':
        [(np.random.rand(10)*10).tolist() for i in range(num_tests)],
    }
    for t in range(num_tests):
        p = get_parameters(params, t)
        numeric_assert("sph_besselj",
                       "sph_besselj",
                       *p,
                       nargout=1,
                       namespace='ars')

def test_sph_bessely():
    num_tests = 10
    params = {
        'n':
        [np.random.randint(10) for i in range(num_tests)],
        'x':
        [(np.random.rand(10)*10).tolist() for i in range(num_tests)],
    }
    for t in range(num_tests):
        p = get_parameters(params, t)
        numeric_assert("sph_bessely",
                       "sph_bessely",
                       *p,
                       nargout=1,
                       namespace='ars')

def test_sph_hankel1():
    num_tests = 10
    params = {
        'n':
        [np.random.randint(10) for i in range(num_tests)],
        'x':
        [(np.random.rand(10)*10).tolist() for i in range(num_tests)],
    }
    for t in range(num_tests):
        p = get_parameters(params, t)
        numeric_assert("sph_hankel1",
                       "sph_hankel1",
                       *p,
                       nargout=1,
                       namespace='ars')

def test_sph_hankel2():
    num_tests = 10
    params = {
        'n':
        [np.random.randint(10) for i in range(num_tests)],
        'x':
        [(np.random.rand(10)*10).tolist() for i in range(num_tests)],
    }
    for t in range(num_tests):
        p = get_parameters(params, t)
        numeric_assert("sph_hankel2",
                       "sph_hankel2",
                       *p,
                       nargout=1,
                       namespace='ars')

def test_dsph_besselj():
    num_tests = 10
    params = {
        'n':
        [np.random.randint(10) for i in range(num_tests)],
        'x':
        [(np.random.rand(10)*10).tolist() for i in range(num_tests)],
    }
    for t in range(num_tests):
        p = get_parameters(params, t)
        numeric_assert("dsph_besselj",
                       "dsph_besselj",
                       *p,
                       nargout=1,
                       namespace='ars')

def test_dsph_bessely():
    num_tests = 10
    params = {
        'n':
        [np.random.randint(10) for i in range(num_tests)],
        'x':
        [(np.random.rand(10)*10).tolist() for i in range(num_tests)],
    }
    for t in range(num_tests):
        p = get_parameters(params, t)
        numeric_assert("dsph_bessely",
                       "dsph_bessely",
                       *p,
                       nargout=1,
                       namespace='ars')

def test_dsph_hankel1():
    num_tests = 10
    params = {
        'n':
        [np.random.randint(10) for i in range(num_tests)],
        'x':
        [(np.random.rand(10)*10).tolist() for i in range(num_tests)],
    }
    for t in range(num_tests):
        p = get_parameters(params, t)
        numeric_assert("dsph_hankel1",
                       "dsph_hankel1",
                       *p,
                       nargout=1,
                       namespace='ars')

def test_dsph_hankel2():
    num_tests = 10
    params = {
        'n':
        [np.random.randint(10) for i in range(num_tests)],
        'x':
        [(np.random.rand(10)*10).tolist() for i in range(num_tests)],
    }
    for t in range(num_tests):
        p = get_parameters(params, t)
        numeric_assert("dsph_hankel2",
                       "dsph_hankel2",
                       *p,
                       nargout=1,
                       namespace='ars')