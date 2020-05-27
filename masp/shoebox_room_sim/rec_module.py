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
#   @file   rec_module.py
#   @author Andrés Pérez-López
#   @date   30/07/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
import copy

import numpy as np
from masp.utils import cart2sph, get_sh, C
from masp.validate_data_types import _validate_echogram_array, _validate_ndarray_2D, _validate_int, _validate_ndarray_1D


def rec_module_mic(echograms, mic_specs):
    """
    Apply microphone directivity gains to a set of given echograms.

    Parameters
    ----------
    echograms : ndarray, dtype = Echogram
        Target echograms. Dimension = (nSrc, nRec)
    mic_specs : ndarray
        Microphone directions and directivity factor. Dimension = (nRec, 4)

    Returns
    -------
    rec_echograms : ndarray, dtype = Echogram
        Echograms subjected to microphone gains. Dimension = (nSrc, nRec)

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    Each row of `mic_specs` is expected to be described as [x, y, z, alpha],
    with (x, y, z) begin the unit vector of the mic orientation.
    `alpha` must be contained in the range [0(dipole), 1(omni)],
    so that directivity is expressed as: d(theta) = a + (1-a)*cos(theta).

    """

    nSrc = echograms.shape[0]
    nRec = echograms.shape[1]
    _validate_echogram_array(echograms)
    _validate_ndarray_2D('mic_specs', mic_specs, shape0=nRec, shape1=C+1)

    mic_vecs = mic_specs[:,:C]
    mic_coeffs = mic_specs[:,-1]

    rec_echograms = copy.copy(echograms)
    # Do nothing if all orders are zeros(omnis)
    if not np.all(mic_coeffs == 1):
        for ns in range(nSrc):
            for nr in range(nRec):
                nRefl = len(echograms[ns, nr].value)
                # Get vectors from source to receiver
                rec_vecs = echograms[ns, nr].coords
                rec_vecs = rec_vecs / np.sqrt(np.sum(np.power(rec_vecs,2), axis=1))[:,np.newaxis]

                mic_gains = mic_coeffs[nr] + (1 - mic_coeffs[nr]) * np.sum(rec_vecs * mic_vecs[nr,:], axis=1)
                rec_echograms[ns, nr].value = echograms[ns, nr].value * mic_gains[:,np.newaxis]

    _validate_echogram_array(rec_echograms)
    return rec_echograms


def rec_module_sh(echograms, sh_orders):
    """
    Apply spherical harmonic directivity gains to a set of given echograms.

    Parameters
    ----------
    echograms : ndarray, dtype = Echogram
        Target echograms. Dimension = (nSrc, nRec)
    sh_orders : int or ndarray, dtype = int
        Spherical harmonic expansion order. Dimension = 1 or (nRec)

    Returns
    -------
    rec_echograms : ndarray, dtype = Echogram
        Echograms subjected to microphone gains. Dimension = (nSrc, nRec)

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    If `sh_orders` is an integer, the given order will be applied to all receivers.

    """

    nSrc = echograms.shape[0]
    nRec = echograms.shape[1]
    _validate_echogram_array(echograms)
    if isinstance(sh_orders, int):
        _validate_int('sh_orders', sh_orders, positive=True)
        sh_orders = sh_orders * np.ones(nRec)
    else:
        _validate_ndarray_1D('sh_orders', sh_orders, size=nRec, positive=True, dtype=int)

    rec_echograms = copy.copy(echograms)
    # Do nothing if all orders are zeros(omnis)
    if not np.all(sh_orders == 0):
        for ns in range(nSrc):
            for nr in range(nRec):
                # Get vectors from source to receiver
                sph = cart2sph(echograms[ns, nr].coords)
                azi = sph[:, 0]
                polar = np.pi / 2 - sph[:, 1]

                sh_gains = get_sh(int(sh_orders[nr]), np.asarray([azi, polar]).transpose(), 'real')
                print(sh_gains)
                # rec_echograms[ns, nr].value = sh_gains * echograms[ns, nr].value[:,np.newaxis]
                rec_echograms[ns, nr].value = sh_gains * echograms[ns, nr].value

    _validate_echogram_array(rec_echograms)
    return rec_echograms
