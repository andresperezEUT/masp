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
#   @file   absorption_module.py
#   @author Andrés Pérez-López
#   @date   30/07/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
import copy

import numpy as np
from .echogram import Echogram

from masp.validate_data_types import _validate_echogram, _validate_ndarray_2D, _validate_ndarray_1D
from masp.utils import C

def apply_absorption(echogram, alpha, limits=None):
    """
    Applies per-band wall absorption to a given echogram.

    Parameters
    ----------
    echogram : Echogram
        Target Echogram
    alpha : ndarray
        Wall absorption coefficients per band. Dimension = (nBands, 6)
    limits : ndarray, optional
        Maximum reflection time per band (RT60). Dimension = (nBands)

    Returns
    -------
    abs_echograms : ndarray, dtype = Echogram
        Array with echograms subject to absorption. Dimension = (1, nBands)

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    `nBands` will be determined by the length of `alpha` first dimension.

    `alpha` must have all values in the range [0,1].

    If 'limits' is not specified, no wall absorption is applied.

    """

    # Validate arguments
    _validate_echogram(echogram)
    _validate_ndarray_2D('abs_wall', alpha, shape1=2*C, norm=True)
    nBands = alpha.shape[0]
    if limits is not None:
        _validate_ndarray_1D('limits', limits, size=nBands, positive=True)

    abs_echograms = np.empty(nBands, dtype=Echogram)

    if limits is None:
        for i in range(nBands):
            abs_echograms[i] = copy.copy(echogram)
    else:
        for nb in range(nBands):
            # Find index of last echogram time element smaller than the given limit
            idx_limit = np.arange(len(echogram.time))[echogram.time < limits[nb]][-1]
            # idx_limit = echogram.time[echogram.time < limits[nb]].size
            abs_echograms[nb] = Echogram(value=echogram.value[:idx_limit+1],
                                         time=echogram.time[:idx_limit+1],
                                         order=echogram.order[:idx_limit+1],
                                         coords=echogram.coords[:idx_limit+1])

    for nb in range(nBands):

        # Absorption coefficients for x, y, z walls per frequency
        a_x = alpha[nb, 0:2]
        a_y = alpha[nb, 2:4]
        a_z = alpha[nb, 4:6]
        # Reflection coefficients
        r_x = np.sqrt(1 - a_x)
        r_y = np.sqrt(1 - a_y)
        r_z = np.sqrt(1 - a_z)

        # Split
        i = abs_echograms[nb].order[:, 0]
        j = abs_echograms[nb].order[:, 1]
        k = abs_echograms[nb].order[:, 2]

        i_even = i[np.remainder(i, 2) == 0]
        i_odd = i[np.remainder(i, 2) != 0]
        i_odd_pos = i_odd[i_odd > 0]
        i_odd_neg = i_odd[i_odd < 0]

        j_even = j[np.remainder(j, 2) == 0]
        j_odd = j[np.remainder(j, 2) != 0]
        j_odd_pos = j_odd[j_odd > 0]
        j_odd_neg = j_odd[j_odd < 0]

        k_even = k[np.remainder(k, 2) == 0]
        k_odd = k[np.remainder(k, 2) != 0]
        k_odd_pos = k_odd[k_odd > 0]
        k_odd_neg = k_odd[k_odd < 0]

        # Find total absorption coefficients by calculating the
        # number of hits on every surface, based on the order per dimension
        abs_x = np.zeros(np.size(abs_echograms[nb].time))
        abs_x[np.remainder(i, 2) == 0] = np.power(r_x[0], (np.abs(i_even) / 2.)) * np.power(r_x[1], (np.abs(i_even) / 2.))
        abs_x[(np.remainder(i, 2) != 0) & (i > 0)] = np.power(r_x[0], np.ceil(i_odd_pos / 2.)) * np.power(r_x[1], np.floor(i_odd_pos / 2.))
        abs_x[(np.remainder(i, 2) != 0) & (i < 0)] = np.power(r_x[0], np.floor(np.abs(i_odd_neg) / 2.)) * np.power(r_x[1], np.ceil(np.abs(i_odd_neg) / 2.))

        abs_y = np.zeros(np.size(abs_echograms[nb].time))
        abs_y[np.remainder(j, 2) == 0] = np.power(r_y[0], (np.abs(j_even) / 2.)) * np.power(r_y[1], (np.abs(j_even) / 2.))
        abs_y[(np.remainder(j, 2) != 0) & (j > 0)] = np.power(r_y[0], np.ceil(j_odd_pos / 2.)) * np.power(r_y[1], np.floor(j_odd_pos / 2.))
        abs_y[(np.remainder(j, 2) != 0) & (j < 0)] = np.power(r_y[0], np.floor(np.abs(j_odd_neg) / 2.)) * np.power(r_y[1], np.ceil(np.abs(j_odd_neg) / 2.))

        abs_z = np.zeros(np.size(abs_echograms[nb].time))
        abs_z[np.remainder(k, 2) == 0] = np.power(r_z[0], (np.abs(k_even) / 2.)) * np.power(r_z[1], (np.abs(k_even) / 2.))
        abs_z[(np.remainder(k, 2) != 0) & (k > 0)] = np.power(r_z[0], np.ceil(k_odd_pos / 2.)) * np.power(r_z[1], np.floor(k_odd_pos / 2,))
        abs_z[(np.remainder(k, 2) != 0) & (k < 0)] = np.power(r_z[0], np.floor(np.abs(k_odd_neg) / 2.)) * np.power(r_z[1], np.ceil(np.abs(k_odd_neg) / 2.))

        s_abs_tot = abs_x * abs_y * abs_z
        # Final amplitude of reflection
        abs_echograms[nb].value = (s_abs_tot * abs_echograms[nb].value.transpose()).transpose()

    return abs_echograms[np.newaxis,:]
