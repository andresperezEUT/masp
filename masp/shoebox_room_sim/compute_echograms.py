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
#   @file   compute_echograms.py
#   @author Andrés Pérez-López
#   @date   30/07/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import numpy as np

from masp.utils import C
from masp.validate_data_types import _validate_ndarray_2D, _validate_ndarray_1D, _validate_int
from .echogram import Echogram
from .image_source_method import ims_coreMtx
from .rec_module import rec_module_mic, rec_module_sh
from .absorption_module import apply_absorption

def compute_echograms_array(room, src, rec, abs_wall, limits):
    """
    Compute the echogram response of a microphone array for a given acoustic scenario.

    Parameters
    ----------
    room : ndarray
        Room dimensions in cartesian coordinates. Dimension = (3) [x, y, z].
    src : ndarray
        Source position in cartesian coordinates. Dimension = (nSrc, 3) [[x, y, z]].
    rec : ndarray
        Receiver position in cartesian coordinates. Dimension = (nRec, 3) [[x, y, z]].
    abs_wall : ndarray
        Wall absorption coefficients per band. Dimension = (nBands, 6)
    limits : ndarray
        Maximum echogram computation time per band.  Dimension = (nBands)

    Returns
    -------
    abs_echograms : ndarray, dtype = Echogram
        Array with rendered echograms. Dimension = (nSrc, nRec, nBands)

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    `src` and `rec` positions are specified from the left ground corner
    of the room, using a left-handed coordinate system.
    `room` refers to the wall dimensions.
    Therefore, their values should be positive and smaller than room dimensions.

              _____    _
             |     |   |
             |     |   |
           x ^     |   | l = r[0]
             |     |   |
             |     |   |
             o---->    -
                  y
             |-----|
                w = r[1]

    `abs_wall` must have all values in the range [0,1].
    `nBands` will be determined by the length of `abs_wall` first dimension.

    TODO: expose type as parameter?, validate return
    """

    nRec = rec.shape[0]
    nSrc = src.shape[0]
    nBands = abs_wall.shape[0]

    _validate_ndarray_1D('room', room, size=C, positive=True)
    _validate_ndarray_2D('src', src, shape1=C, positive=True)
    _validate_ndarray_2D('rec', rec, shape1=C, positive=True)
    _validate_ndarray_2D('abs_wall', abs_wall, shape1=2*C, positive=True)
    _validate_ndarray_1D('limits', limits, positive=True, size=nBands)

    # Limit the RIR by reflection order or by time-limit
    type = 'maxTime'

    echograms = np.empty((nSrc, nRec), dtype=Echogram)
    # Compute echogram due to pure propagation (frequency-independent)
    for ns in range(nSrc):
        for nr in range(nRec):
            print('Compute echogram: Source ' + str(ns) + ' - Receiver ' + str(nr))
            # Compute echogram
            echograms[ns, nr] = ims_coreMtx(room, src[ns,:], rec[nr,:], type, np.max(limits))

    abs_echograms = np.empty((nSrc, nRec, nBands), dtype=Echogram)
    # Apply boundary absorption
    for ns in range(nSrc):
        for nr in range(nRec):
            print('Apply absorption: Source ' + str(ns) + ' - Receiver ' + str(nr))
            # Compute echogram
            abs_echograms[ns, nr] = apply_absorption(echograms[ns, nr], abs_wall, limits)

    # return abs_echograms, echograms
    return abs_echograms



def compute_echograms_mic(room, src, rec, abs_wall, limits, mic_specs):
    """
    Compute the echogram response of individual microphones for a given acoustic scenario.

    Parameters
    ----------
    room : ndarray
        Room dimensions in cartesian coordinates. Dimension = (3) [x, y, z].
    src : ndarray
        Source position in cartesian coordinates. Dimension = (nSrc, 3) [[x, y, z]].
    rec : ndarray
        Receiver position in cartesian coordinates. Dimension = (nRec, 3) [[x, y, z]].
    abs_wall : ndarray
        Wall absorption coefficients per band. Dimension = (nBands, 6)
    limits : ndarray
        Maximum echogram computation time per band.  Dimension = (nBands)
    mic_specs : ndarray
        Microphone directions and directivity factor. Dimension = (nRec, 4)

    Returns
    -------
    abs_echograms : ndarray, dtype = Echogram
        Array with rendered echograms. Dimension = (nSrc, nRec, nBands)

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    `src` and `rec` positions are specified from the left ground corner
    of the room, using a left-handed coordinate system.
    `room` refers to the wall dimensions.
    Therefore, their values should be positive and smaller than room dimensions.

              _____    _
             |     |   |
             |     |   |
           x ^     |   | l = r[0]
             |     |   |
             |     |   |
             o---->    -
                  y
             |-----|
                w = r[1]

    `abs_wall` must have all values in the range [0,1].
    `nBands` will be determined by the length of `abs_wall` first dimension.

    Each row of `mic_specs` is expected to be described as [x, y, z, alpha],
    with (x, y, z) begin the unit vector of the mic orientation.
    `alpha` must be contained in the range [0(dipole), 1(omni)],
    so that directivity is expressed as: d(theta) = a + (1-a)*cos(theta).


    TODO: expose type as parameter?, validate return
    """

    nRec = rec.shape[0]
    nSrc = src.shape[0]
    nBands = abs_wall.shape[0]

    _validate_ndarray_1D('room', room, size=C, positive=True)
    _validate_ndarray_2D('src', src, shape1=C, positive=True)
    _validate_ndarray_2D('rec', rec, shape1=C, positive=True)
    _validate_ndarray_2D('abs_wall', abs_wall, shape1=2*C, positive=True)
    _validate_ndarray_1D('limits', limits, positive=True, size=nBands)
    _validate_ndarray_2D('mic_specs', mic_specs, shape0=nRec, shape1=C+1)

    # Limit the RIR by reflection order or by time-limit
    type = 'maxTime'

    # Compute echogram due to pure propagation (frequency-independent)
    echograms = np.empty((nSrc, nRec), dtype=Echogram)
    for ns in range(nSrc):
        for nr in range(nRec):
            print('Compute echogram: Source ' + str(ns) + ' - Receiver ' + str(nr))
            # Compute echogram
            echograms[ns, nr] = ims_coreMtx(room, src[ns,:], rec[nr,:], type, np.max(limits))

    print('Apply receiver direcitivites')
    rec_echograms = rec_module_mic(echograms, mic_specs)

    abs_echograms = np.empty((nSrc, nRec, nBands), dtype=Echogram)
    # Apply boundary absorption
    for ns in range(nSrc):
        for nr in range(nRec):
            print('Apply absorption: Source ' + str(ns) + ' - Receiver ' + str(nr))
            # Compute echogram
            abs_echograms[ns, nr] = apply_absorption(rec_echograms[ns, nr], abs_wall, limits)

    # return abs_echograms, rec_echograms, echograms
    return abs_echograms


def compute_echograms_sh(room, src, rec, abs_wall, limits, sh_orders):
    """
    Compute the echogram response of individual microphones for a given acoustic scenario,
    in the spherical harmonic domain.

    Parameters
    ----------
    room : ndarray
        Room dimensions in cartesian coordinates. Dimension = (3) [x, y, z].
    src : ndarray
        Source position in cartesian coordinates. Dimension = (nSrc, 3) [[x, y, z]].
    rec : ndarray
        Receiver position in cartesian coordinates. Dimension = (nRec, 3) [[x, y, z]].
    abs_wall : ndarray
        Wall absorption coefficients per band. Dimension = (nBands, 6)
    limits : ndarray
        Maximum echogram computation time per band.  Dimension = (nBands)
    sh_orders : int or ndarray, dtype = int
        Spherical harmonic expansion order. Dimension = 1 or (nRec)

    Returns
    -------
    abs_echograms : ndarray, dtype = Echogram
        Array with rendered echograms. Dimension = (nSrc, nRec, nBands)

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    `src` and `rec` positions are specified from the left ground corner
    of the room, using a left-handed coordinate system.
    `room` refers to the wall dimensions.
    Therefore, their values should be positive and smaller than room dimensions.

              _____    _
             |     |   |
             |     |   |
           x ^     |   | l = r[0]
             |     |   |
             |     |   |
             o---->    -
                  y
             |-----|
                w = r[1]

    `abs_wall` must have all values in the range [0,1].
    `nBands` will be determined by the length of `abs_wall` first dimension.

    If `sh_orders` is an integer, the given order will be applied to all receivers.
    'nRec' will be determined by the length of `rec` first dimension.

    TODO: expose type as parameter?, validate return
    """

    nRec = rec.shape[0]
    nSrc = src.shape[0]
    nBands = abs_wall.shape[0]

    _validate_ndarray_1D('room', room, size=C, positive=True)
    _validate_ndarray_2D('src', src, shape1=C, positive=True)
    _validate_ndarray_2D('rec', rec, shape1=C, positive=True)
    _validate_ndarray_2D('abs_wall', abs_wall, shape1=2*C, positive=True)
    _validate_ndarray_1D('limits', limits, positive=True, size=nBands)
    if isinstance(sh_orders, int):
        _validate_int('sh_orders', sh_orders, positive=True)
        sh_orders = sh_orders * np.ones(nRec, dtype=int)
    else:
        _validate_ndarray_1D('sh_orders', sh_orders, size=nRec, positive=True, dtype=int)

    # Limit the RIR by reflection order or by time-limit
    type = 'maxTime'
    # Compute echogram due to pure propagation (frequency-independent)
    echograms = np.empty((nSrc, nRec), dtype=Echogram)
    for ns in range(nSrc):
        for nr in range(nRec):
            print('Compute echogram: Source ' + str(ns) + ' - Receiver ' + str(nr))
            # Compute echogram
            echograms[ns, nr] = ims_coreMtx(room, src[ns,:], rec[nr,:], type, np.max(limits))

    print('Apply SH directivites')
    rec_echograms = rec_module_sh(echograms, sh_orders)

    abs_echograms = np.empty((nSrc, nRec, nBands), dtype=Echogram)
    # Apply boundary absorption
    for ns in range(nSrc):
        for nr in range(nRec):
            print ('Apply absorption: Source ' + str(ns) + ' - Receiver ' + str(nr))
            # Compute echogram
            abs_echograms[ns, nr] = apply_absorption(rec_echograms[ns, nr], abs_wall, limits)

    # return abs_echograms, rec_echograms, echograms
    return abs_echograms

