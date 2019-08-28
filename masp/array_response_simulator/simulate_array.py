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
#   @file   simulate_array.py
#   @author Andrés Pérez-López
#   @date   22/08/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import numpy as np
import scipy.special
import masp.array_response_simulator as asr
import masp.utils
from masp.validate_data_types import _validate_int, _validate_ndarray_2D, _validate_string, _validate_float, \
    _validate_ndarray_1D


def simulate_sph_array(N_filt, mic_dirs_rad, src_dirs_rad, arrayType, R, N_order, fs, dirCoef=None):
    """
    Simulate the impulse responses of a spherical array.

    Parameters
    ----------
    N_filt : int
        Number of frequencies where to compute the response. It must be even.
    mic_dirs_rad: ndarray
        Directions of microphone capsules, in radians.
        Expressed in [azi, ele] pairs. Dimension = (N_mic, C-1).
    src_dirs_rad: ndarray
        Direction of arrival of the indicent plane waves, in radians.
         Expressed in [azi, ele] pairs. Dimension = (N_doa, C-1).
    arrayType: str
        'open', 'rigid' or 'directional'.
        Target sampling rate
    R: float
        Radius of the array sphere, in meter.
    N_order: int
        Maximum spherical harmonic expansion order.
    fs: int
        Sample rate.
    dirCoef: float, optional
        Directivity coefficient of the sensors. Default to None.

    Returns
    -------
    h_mic: ndarray
        Computed IRs in time-domain. Dimension = (N_filt, N_mic, N_doa).
    H_mic: ndarray, dtype='complex'
        Frequency responses of the computed IRs. Dimension = (N_filt//2+1, N_mic, N_doa).

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    This method computes the impulse responses of the microphones of a
    spherical microphone array for the given directions of incident plane waves.
    The array type can be either 'open' for omnidirectional microphones in
    an open setup, 'rigid' for omnidirectional microphones mounted on a sphere, or
    'directional' for an array of first-order directional microphones determined
    by `dirCoef`.
    """

    _validate_int('N_filt', N_filt, positive=True, parity='even')
    _validate_ndarray_2D('mic_dirs_rad', mic_dirs_rad, shape1=masp.utils.C-1)
    _validate_ndarray_2D('src_dirs_rad', src_dirs_rad, shape1=masp.utils.C-1)
    _validate_string('arrayType', arrayType, choices=['open', 'rigid', 'directional'])
    _validate_float('R', R, positive=True)
    _validate_int('N_order', N_order, positive=True)
    _validate_int('fs', fs, positive=True)
    if arrayType is 'directional':
        if dirCoef is None:
            raise ValueError('dirCoef must be defined in the directional case.')
        _validate_float('dirCoef', dirCoef)

    # Compute the frequency-dependent part of the microphone responses (radial dependence)
    f = np.arange(N_filt//2+1) * fs / N_filt
    c = 343.
    kR = 2*np.pi*f*R/c
    b_N = asr.sph_modal_coefs(N_order, kR, arrayType, dirCoef)

    # Handle Nyquist for real impulse response
    temp = b_N.copy()
    temp[-1,:] = np.real(temp[-1,:])
    # Create the symmetric conjugate negative frequency response for a real time-domain signal
    b_Nt = np.real(np.fft.fftshift(np.fft.ifft(np.append(temp, np.conj(temp[-2:0:-1,:]), axis=0), axis=0), axes=0))

    # Compute angular-dependent part of the microphone responses
    # Unit vectors of DOAs and microphones
    N_doa = src_dirs_rad.shape[0]
    N_mic = mic_dirs_rad.shape[0]
    U_doa = masp.utils.sph2cart(src_dirs_rad[:, 0], src_dirs_rad[:, 1], 1)
    U_mic = masp.utils.sph2cart(mic_dirs_rad[:, 0], mic_dirs_rad[:, 1], 1)

    h_mic = np.zeros((N_filt, N_mic, N_doa))
    H_mic = np.zeros((N_filt // 2 + 1, N_mic, N_doa), dtype='complex')

    for i in range(N_doa):
        cosangle = np.dot(U_mic, U_doa[i,:])
        P = np.zeros((N_order + 1, N_mic))
        for n in range(N_order+1):
            for mic in range(N_mic):
                # The Legendre polynomial gives the angular dependency
                Pn = scipy.special.lpmn(n,n,cosangle[mic])[0][0,-1]
                P[n , mic] = (2 * n + 1) / (4 * np.pi) * Pn
        h_mic[:, :, i] = np.matmul(b_Nt, P)
        H_mic[:, :, i] = np.matmul(b_N, P)

    return h_mic, H_mic


def simulate_cyl_array(N_filt, mic_dirs_rad, src_dirs_rad, arrayType, R, N_order, fs):
    """
    Simulate the impulse responses of a cylindrical array.

    Parameters
    ----------
    N_filt : int
        Number of frequencies where to compute the response. It must be even.
    mic_dirs_rad: ndarray
        Directions of microphone capsules, in radians. Dimension = (N_mic).
    src_dirs_rad: ndarray
        Direction of arrival of the indicent plane waves, in radians. Dimension = (N_doa).
    arrayType: str
        'open' or 'rigid'.
        Target sampling rate
    R: float
        Radius of the array cylinder, in meter.
    N_order: int
        Maximum cylindrical harmonic expansion order.
    fs: int
        Sample rate.

    Returns
    -------
    h_mic: ndarray
        Computed IRs in time-domain. Dimension = (N_filt, N_mic, N_doa).
    H_mic: ndarray, dtype='complex'
        Frequency responses of the computed IRs. Dimension = (N_filt//2+1, N_mic, N_doa).

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    This method computes the impulse responses of the microphones of a
    cylindrical microphone array for the given directions of incident plane waves.
    The array type can be either 'open' for omnidirectional microphones in
    an open setup, or 'rigid' for omnidirectional microphones mounted on a cylinder.

    """

    _validate_int('N_filt', N_filt, positive=True, parity='even')
    _validate_ndarray_1D('mic_dirs_rad', mic_dirs_rad)
    _validate_ndarray_1D('src_dirs_rad', src_dirs_rad)
    _validate_string('arrayType', arrayType, choices=['open', 'rigid'])
    _validate_float('R', R, positive=True)
    _validate_int('N_order', N_order, positive=True)
    _validate_int('fs', fs, positive=True)

    # Compute the frequency-dependent part of the microphone responses (radial dependence)
    f = np.arange(N_filt//2+1) * fs / N_filt
    c = 343.
    kR = 2*np.pi*f*R/c
    b_N = asr.cyl_modal_coefs(N_order, kR, arrayType)

    # Handle Nyquist for real impulse response
    temp = b_N.copy()
    temp[-1,:] = np.real(temp[-1,:])
    # Create the symmetric conjugate negative frequency response for a real time-domain signal
    b_Nt = np.real(np.fft.fftshift(np.fft.ifft(np.append(temp, np.conj(temp[-2:0:-1,:]), axis=0), axis=0), axes=0))

    # Compute angular-dependent part of the microphone responses
    # Unit vectors of DOAs and microphones
    N_doa = src_dirs_rad.shape[0]
    N_mic = mic_dirs_rad.shape[0]
    h_mic = np.zeros((N_filt, N_mic, N_doa))
    H_mic = np.zeros((N_filt // 2 + 1, N_mic, N_doa), dtype='complex')

    for i in range(N_doa):
        angle = mic_dirs_rad - src_dirs_rad[i]
        C = np.zeros((N_order + 1, N_mic))
        for n in range(N_order+1):
            # Jacobi-Anger expansion
            if n == 0:
                C[n, :] = np.ones(angle.shape)
            else:
                C[n, :] = 2 * np.cos(n*angle)
        h_mic[:, :, i] = np.matmul(b_Nt, C)
        H_mic[:, :, i] = np.matmul(b_N, C)

    return h_mic, H_mic


