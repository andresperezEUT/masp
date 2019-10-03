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
#   @file   array_sht_filters.py
#   @author Andrés Pérez-López
#   @date   02/10/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import numpy as np
import warnings
from masp.utils import c
from masp.array_response_simulator import sph_modal_coefs
from masp.validate_data_types import _validate_int, _validate_float

def array_sht_filters_theory_radInverse(R, nMic, order_sht, Lfilt, fs, amp_threshold):
    """
    Generate SHT filters based on theoretical responses (regularized radial inversion)

    Parameters
    ----------
    R : float
        Microphone array radius, in meter.
    nMic : int
        Number of microphone capsules.
    order_sht : int
        Spherical harmonic transform order.
    Lfilt : int
        Number of FFT points for the output filters. It must be even.
    fs : int
        Sample rate for the output filters.
    amp_threshold : float
         Max allowed amplification for filters, in dB.

    Returns
    -------
    h_filt : ndarray
        Impulse responses of the filters. Dimension = (Lfilt, order_sht+1).
    H_filt: ndarray
        Frequency-domain filters. Dimension = (Lfilt//2+1, order_sht+1).

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.
    UserWarning: if `nMic` not big enough for the required sht order.

    Notes
    -----
    Generate the filters to convert microphone signals from a spherical
    microphone array to SH signals, based on an ideal theoretical model of
    the array. The filters are generated from an inversion of the radial
    components of the response, neglecting spatial aliasing effects and
    non-ideal arrangements of the microphones. One filter is shared by all
    SH signals of the same order in this case.
    Here this single channel inversion problem is done with a constraint on
    max allowed amplification regularization, using Tikhonov
    regularization, similar e.g. to:

        Moreau, S., Daniel, J., Bertet, S., 2006,
        3D sound field recording with higher order ambisonics-objective
        measurements and validation of spherical microphone.
        In Audio Engineering Society Convention 120.

    """

    _validate_float('R', R, positive=True)
    _validate_int('nMic', nMic, positive=True)
    _validate_int('order_sht', order_sht, positive=True)
    _validate_int('Lfilt', Lfilt, positive=True, parity='even')
    _validate_int('fs', fs, positive=True)
    _validate_float('amp_threshold', amp_threshold)

    f = np.arange(Lfilt // 2 + 1) * fs / Lfilt

    # Adequate sht order to the number of microphones
    if order_sht > np.sqrt(nMic) - 1:
        order_sht = int(np.floor(np.sqrt(nMic) - 1))
        warnings.warn("Set order too high for the number of microphones, should be N<=np.sqrt(Q)-1. Auto set to "+str(order_sht), UserWarning)

    # Modal responses
    kR = 2 * np.pi * f * R / c
    bN = sph_modal_coefs(order_sht, kR, 'rigid') / (4 * np.pi)

    # Encoding matrix with regularization
    a_dB = amp_threshold
    alpha = complex(np.sqrt(nMic) * np.power(10, a_dB / 20))  # Explicit casting to allow negative sqrt (a_dB < 0)
    beta = np.sqrt((1 - np.sqrt(1 - 1 / np.power(alpha,2))) / (1 + np.sqrt(1 - 1 / np.power(alpha,2))))  # Moreau & Daniel
    # Regularized single channel equalization filters per order
    H_filt = (bN).conj() / (np.power(np.abs(bN), 2) + np.power(beta, 2) * np.ones((Lfilt // 2 + 1, order_sht + 1)))

    # Time domain filters
    temp = H_filt.copy()
    temp[-1,:] = np.real(temp[-1,:])
    # Create the symmetric conjugate negative frequency response for a real time-domain signal
    h_filt = np.real(np.fft.fftshift(np.fft.ifft(np.append(temp, np.conj(temp[-2:0:-1,:]), axis=0), axis=0), axes=0))

    return h_filt, H_filt
    # TODO: check return order
    # return H_filt, h_filt
