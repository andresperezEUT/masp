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
#   @file   sph_modal_coefs.py
#   @author Andrés Pérez-López
#   @date   22/08/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import numpy as np
from masp.array_response_simulator import sph_besselj


def sph_modal_coefs(N, kr, arrayType, dirCoeff):
    """
    %SPHMODALCOEFFS Modal coefficients for rigid or open spherical array
    %
    %   N: maximum order
    %   kr: wavenumber-radius product
    %   arrayType: {'open','rigid','directional'} open for open array of
    %               omnidirectional sensors, rigid for sensors mounted on a
    %               rigid baffle, directional for an array of first-order
    %               directional microphones determined by the dirCoeff
    %   dirCoeff:   relevant only in the 'directional' type, dirCoeff ranges
    %               from 0 (omni) to 1 (dipole), where for example 0.5 is a
    %               cardioid sensor. In the 0 case it is equivalent to an open
    %               array of omnis. The first order directivity function is
    %               defined as d(theta) = dirCoeff + (1-dirCoeff)*cos(theta)

    :param N: INT
    :param kr: POSITIVE
    :param arrayType:
    :param dirCoeff:
    :return:

    TODO
    """

    b_N = np.zeros((kr.size, N+1))

    for n in range(N):

        if arrayType is 'open':
            b_N[:,n+1] = 4 * np.pi * np.power(1j,n) * sph_besselj(n, kr)

        elif arrayType is 'rigid':
            jn = sph_besselj(n, kr)
            jnprime = dsph_besselj(n, kr);


