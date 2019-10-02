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
#   @file   test_utils.py
#   @author Andrés Pérez-López
#   @date   12/08/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
from masp import cart2sph, sph2cart, elev2incl, incl2elev
from masp.tests.convenience_test_methods import *
import pytest



def test_cart2sph():

    # Regular 2D case
    cart = np.asarray([
        [1., 0., 0.],
        [-1., 0., 0.],
        [0., 1., 0.],
        [0., -1., 0.],
        [0., 0., 1.],
        [0., 0., -1.],
        [1., 1., 1.],
        [-1., -1., -1.]
    ])

    sph = np.asarray([
        [0., 0., 1],
        [np.pi, 0., 1.],
        [np.pi/2, 0., 1.],
        [-np.pi/2, 0., 1.],
        [0., np.pi/2, 1.],
        [0., -np.pi/2, 1.],
        [np.pi/4, np.arcsin(1. / np.sqrt(3)), np.sqrt(3)],
        [-3*np.pi/4, -np.arcsin(1. / np.sqrt(3)), np.sqrt(3)],
    ])

    assert np.allclose(sph, cart2sph(cart))

    # Specific 1D case
    cart = np.asarray([1., 0., 0.])
    sph = np.asarray([0., 0., 1])
    assert cart2sph(cart).ndim == 1
    assert np.allclose(sph, cart2sph(cart))

    # Error with higher dimensions
    with pytest.raises(ValueError):
        cart2sph(np.ones((1, 1, 1)))


def test_sph2cart():

    sph = np.asarray([
        [0., 0., 1],
        [np.pi, 0., 1.],
        [np.pi/2, 0., 1.],
        [-np.pi/2, 0., 1.],
        [0., np.pi/2, 1.],
        [0., -np.pi/2, 1.],
        [np.pi/4, np.arcsin(1. / np.sqrt(3)), np.sqrt(3)],
        [-3*np.pi/4, -np.arcsin(1. / np.sqrt(3)), np.sqrt(3)],
    ])

    cart = np.asarray([
        [1., 0., 0.],
        [-1., 0., 0.],
        [0., 1., 0.],
        [0., -1., 0.],
        [0., 0., 1.],
        [0., 0., -1.],
        [1., 1., 1.],
        [-1., -1., -1.]
    ])

    assert np.allclose(cart, sph2cart(sph))

    # Specific 1D case
    sph = np.asarray([0., 0., 1])
    cart = np.asarray([1., 0., 0.])
    assert sph2cart(sph).ndim == 1
    assert np.allclose(cart, sph2cart(sph))

    # Error with higher dimensions
    with pytest.raises(ValueError):
        sph2cart(np.ones((1, 1, 1)))


def test_elev2incl():

    # shape1=2
    elev = np.asarray([
        [np.pi/2, np.pi/2],
        [np.pi/4, np.pi/4],
        [0., 0.],
        [-np.pi/4, -np.pi/4],
        [-np.pi/2, -np.pi/2],
    ])
    incl = np.asarray([
        [np.pi/2, 0],
        [np.pi/4, np.pi/4],
        [0., np.pi/2],
        [-np.pi/4, 3*np.pi/4],
        [-np.pi/2, np.pi],
    ])
    assert np.allclose(incl, elev2incl(elev))

    # shape1=3
    elev = np.asarray([
        [np.pi/2, np.pi/2, 1.],
        [np.pi/4, np.pi/4, 1.],
        [0., 0., 1.],
        [-np.pi/4, -np.pi/4, 1.],
        [-np.pi/2, -np.pi/2, 1.],
    ])
    incl = np.asarray([
        [np.pi/2, 0, 1.],
        [np.pi/4, np.pi/4, 1.],
        [0., np.pi/2, 1.],
        [-np.pi/4, 3*np.pi/4, 1.],
        [-np.pi/2, np.pi, 1.],
    ])
    assert np.allclose(incl, elev2incl(elev))

    # Error with other shapes
    wrong_values = [
        np.ones((1, 1, 1)),
        np.ones(3),
        np.zeros((5, 1)),
        np.zeros((3, 4)),
    ]
    for wv in wrong_values:
        print(wv)
        with pytest.raises(ValueError):
            elev2incl(wv)

def test_incl2elev():
    # shape1=2
    incl = np.asarray([
        [np.pi/2, 0],
        [np.pi/4, np.pi/4],
        [0., np.pi/2],
        [-np.pi/4, 3*np.pi/4],
        [-np.pi/2, np.pi],
    ])
    elev = np.asarray([
        [np.pi/2, np.pi/2],
        [np.pi/4, np.pi/4],
        [0., 0.],
        [-np.pi/4, -np.pi/4],
        [-np.pi/2, -np.pi/2],

    ])
    assert np.allclose(elev, incl2elev(incl))

    # shape1=3
    incl = np.asarray([
        [np.pi/2, 0, 1.],
        [np.pi/4, np.pi/4, 1.],
        [0., np.pi/2, 1.],
        [-np.pi/4, 3*np.pi/4, 1.],
        [-np.pi/2, np.pi, 1.],
    ])
    elev = np.asarray([
        [np.pi/2, np.pi/2, 1.],
        [np.pi/4, np.pi/4, 1.],
        [0., 0., 1.],
        [-np.pi/4, -np.pi/4, 1.],
        [-np.pi/2, -np.pi/2, 1.],
    ])
    assert np.allclose(elev, incl2elev(incl))

    # Error with other shapes
    wrong_values = [
        np.ones((1, 1, 1)),
        np.ones(3),
        np.zeros((5, 1)),
        np.zeros((3, 4)),
    ]
    for wv in wrong_values:
        with pytest.raises(ValueError):
            incl2elev(wv)


def test_get_sh():
    num_tests = 10
    params = {
        'N':
        [np.random.randint(1,10) for i in range(num_tests)],
        'dirs':
        [(np.random.rand(np.random.randint(1,10),2)*[2*np.pi, np.pi]).tolist() for i in range(num_tests)],
        'basisType':
        # [random.choice(['real', 'complex']) for i in range(num_tests)], # TODO: to be uncommented
        ['real' for i in range(num_tests)],
    }
    for t in range(num_tests):
        p = get_parameters(params, t)
        numeric_assert("getSH",
                       "get_sh",
                       *p,
                       nargout=1)

def test_lagrange():
    num_tests = 10
    ns = [np.random.randint(1,101) for i in range(num_tests)]
    params = {
        'N':
        ns,
        'delays':
        [(np.linspace(0, 1, np.random.randint(3, 11)) + (ns[i] / 2)).tolist() for i in range(num_tests)]        # [(np.random.rand(np.random.randint(1,10),2)*[2*np.pi, np.pi]).tolist() for i in range(num_tests)],
    }
    for t in range(num_tests):
        p = get_parameters(params, t)
        numeric_assert("lagrange",
                       "lagrange",
                       *p,
                       nargout=1)


def test_isLambda():

    wrong_values = [1, '1', True, 2.3, 1e4, 3j, [1], None, np.nan, np.inf, np.asarray([0.5])]
    for wv in wrong_values:
        assert not masp.isLambda(wv)

    lambda_expr = lambda a : a + 10
    assert masp.isLambda(lambda_expr)


def test_check_cond_number():
    num_tests = 10
    params = {
        'N':
        [np.random.randint(1,10) for i in range(num_tests)],
        'dirs':
        [(np.random.rand(np.random.randint(1,10),C-1)*[2*np.pi, np.pi]).tolist() for i in range(num_tests)],
        'basisType':
        # [random.choice(['real', 'complex']) for i in range(num_tests)], # TODO: to be uncommented
        ['real' for i in range(num_tests)],
        'W':
        [None for i in range(num_tests)]  # TODO not implemented
    }
    for t in range(num_tests):
        p = get_parameters(params, t)
        numeric_assert("checkCondNumberSHT",
                       "check_cond_number_sht",
                       *p,
                       nargout=1,
                       rtol=1e2)  # Big rtol due to numeric errors in matrix condition computation
