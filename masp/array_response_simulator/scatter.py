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
#   @file   scatter.py
#   @author Andrés Pérez-López
#   @date   29/08/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from scipy.special import jv, jvp, hankel2, h2vp, spherical_jn, lpmn
import numpy as np
import masp.utils
import masp.array_response_simulator as asr
from masp.validate_data_types import _validate_ndarray_2D, _validate_float, _validate_int, _validate_ndarray_1D


def spherical_scatterer(mic_dirs_rad, src_dirs_rad, R, N_order, N_filt, fs):
    """
    Compute the pressure due to a spherical scatterer

    The function computes the impulse responses of the pressure measured
    at some points in the field with a spherical rigid scatterer centered
    at the origin and due to incident plane waves.

    Parameters
    ----------
    mic_dirs_rad: ndarray
        Position of microphone capsules. Dimension = (N_mic, C).
        Positions are expected in radians, expressed in triplets [azimuth, elevation, distance].
    src_dirs_rad: ndarray
        Direction of arrival of the indicent plane waves. Dimension = (N_doa, C-1).
        Directions are expected in radians, expressed in pairs [azimuth, elevation].
    R: float
        Radius of the array sphere, in meter.
    N_order: int
        Maximum cylindrical harmonic expansion order.
    N_filt : int
        Number of frequencies where to compute the response. It must be even.
    fs: int
        Sample rate.

    Returns
    -------
    h_mic: ndarray
        Computed IRs in time-domain. Dimension = (N_filt, Nmic, Ndoa)
    H_mic: ndarray, dtype='complex'
        Frequency responses of the computed IRs. Dimension = (N_filt//2+1, N_mic, N_doa).

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    """

    _validate_ndarray_2D('mic_dirs_rad', mic_dirs_rad, shape1=masp.C)
    _validate_ndarray_2D('src_dirs_rad', src_dirs_rad, shape1=masp.C-1)
    _validate_float('R', R, positive=True)
    _validate_int('N_order', N_order, positive=True)
    _validate_int('N_filt', N_filt, positive=True, parity='even')
    _validate_int('fs', fs, positive=True)
    if np.any(mic_dirs_rad[:,2] < R):
        raise ValueError('mic_dirs_rad: The distance of the measurement point cannot be less than the radius:' + str(R))

    # Compute the frequency-dependent part of the microphone responses (radial dependence)
    K = N_filt // 2 + 1
    f = np.arange(K) * fs / N_filt
    kR = 2*np.pi*f*R/masp.c
    N_mic = mic_dirs_rad.shape[0]
    N_doa = src_dirs_rad.shape[0]

    # Check if all microphones at same radius
    same_radius = np.sum(mic_dirs_rad[1:, 2] - mic_dirs_rad[:-1, 2]) == 0
    if same_radius:
        # Spherical modal coefs for rigid sphere
        b_N = np.zeros((K, N_order + 1), dtype='complex')
        r = mic_dirs_rad[0, 2]
        kr = 2 * np.pi * f * r / masp.c

        # Similar to the sph_modal_coefs for the rigid case
        for n in range(N_order + 1):
            jn = spherical_jn(n, kr)
            jnprime = spherical_jn(n, kR, derivative=True)
            hn = asr.sph_hankel2(n, kr)
            hnprime = asr.dsph_hankel2(n, kR)
            b_N[:, n] = (2 * n + 1) * np.power(1j, n) * (jn - (jnprime / hnprime) * hn)
    else:
        # Spherical modal coefs for rigid sphere, but at different distances
        b_N = np.zeros((K, N_order + 1, N_mic), dtype='complex')
        for nm in range(N_mic):
            r = mic_dirs_rad[nm, 2]
            kr = 2 * np.pi * f * r / masp.c

            # Similar to the sph_modal_coefs for the rigid case
            for n in range(N_order+1):
                jn = spherical_jn(n, kr)
                jnprime = spherical_jn(n, kR, derivative=True)
                hn = asr.sph_hankel2(n, kr)
                hnprime = asr.dsph_hankel2(n, kR)
                b_N[:, n, nm] = (2*n+1) * np.power(1j, n) * (jn - (jnprime / hnprime) * hn)

    # Avoid NaNs for very high orders, instead of (very) very small values
    b_N[np.isnan(b_N)] = 0.

    # Compute angular-dependent part of the microphone responses
    H_mic = np.zeros((K, N_mic, N_doa), dtype='complex')
    for nd in range(N_doa):
        # Unit vectors of DOAs and microphones
        azi0 = src_dirs_rad[nd, 0]
        elev0 = src_dirs_rad[nd, 1]
        azi = mic_dirs_rad[:, 0]
        elev = mic_dirs_rad[:, 1]
        cosAlpha = np.sin(elev)*np.sin(elev0) + np.cos(elev)*np.cos(elev0)*np.cos(azi-azi0)

        P_N = np.zeros((N_order+1, N_mic))
        for n in range(N_order+1):
            for nm in range(N_mic):
                P_N[n, nm] = lpmn(n, n, cosAlpha[nm])[0][0, -1]

        # Accumulate across orders
        if same_radius:
            H_mic[:, :, nd] = np.matmul(b_N, P_N)
        else:
            for nm in range(N_mic):
                H_mic[:, nm, nd] = np.matmul(b_N[:, :, nm], P_N[:, nm])

    # Handle Nyquist for real impulse response
    tempH_mic = H_mic.copy()
    # TODO: in `simulate_sph_array()` it was real, not abs. Why?
    tempH_mic[-1, :] = np.abs(tempH_mic[-1, :])
    # Conjugate ifft and fftshift for causal IR
    h_mic = np.real(np.fft.fftshift(np.fft.ifft(np.append(tempH_mic, np.conj(tempH_mic[-2:0:-1,:]), axis=0), axis=0), axes=0))

    return h_mic, H_mic


def cylindrical_scatterer(mic_dirs_rad, src_dirs_rad, R, N_order, N_filt, fs):
    """
    Compute the pressure due to a cylindrical scatterer

    The function computes the impulse responses of the pressure measured
    at some points in the field with a cylindrical rigid scatterer centered
    at the origin and due to incident plane waves.

    Parameters
    ----------
    mic_dirs_rad: ndarray
        Position of microphone capsules. Dimension = (N_mic, C-1).
        Positions are expected in radians, expressed in pairs [azimuth, distance].
    src_dirs_rad: ndarray
        Direction of arrival of the indicent plane waves. Dimension = (N_doa).
        Directions (azimuths) are expected in radians.
    R: float
        Radius of the array sphere, in meter.
    N_order: int
        Maximum cylindrical harmonic expansion order.
    N_filt : int
        Number of frequencies where to compute the response. It must be even.
    fs: int
        Sample rate.

    Returns
    -------
    h_mic: ndarray
        Computed IRs in time-domain. Dimension = (N_filt, Nmic, Ndoa)
    H_mic: ndarray, dtype='complex'
        Frequency responses of the computed IRs. Dimension = (N_filt//2+1, N_mic, N_doa).

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    """

    _validate_ndarray_2D('mic_dirs_rad', mic_dirs_rad, shape1=masp.C-1)
    _validate_ndarray_1D('src_dirs_rad', src_dirs_rad)
    _validate_float('R', R, positive=True)
    _validate_int('N_order', N_order, positive=True)
    _validate_int('N_filt', N_filt, positive=True, parity='even')
    _validate_int('fs', fs, positive=True)
    if np.any(mic_dirs_rad[:, 1] < R):
        raise ValueError('mic_dirs_rad: The distance of the measurement point cannot be less than the radius:' + str(R))

    K = N_filt // 2 + 1
    f = np.arange(K) * fs / N_filt
    kR = 2 * np.pi * f * R / masp.c
    N_mic = mic_dirs_rad.shape[0]
    N_doa = src_dirs_rad.size

    # Check if all microphones at same radius
    same_radius = np.sum(mic_dirs_rad[1:, 1] - mic_dirs_rad[:-1, 1]) == 0
    if same_radius:
        # Cylindrical modal coefs for rigid sphere
        b_N = np.zeros((K, N_order + 1), dtype='complex')
        r = mic_dirs_rad[0, 1]
        kr = 2 * np.pi * f * r / masp.c

        # Similar to the cyl_modal_coefs for the rigid case
        for n in range(N_order + 1):
            jn = jv(n, kr)
            jnprime = jvp(n, kR, 1)
            hn = hankel2(n, kr)
            hnprime = h2vp(n, kR, 1)
            b_N[:, n] = np.power(1j, n) * (jn - (jnprime / hnprime) * hn)

    else:
        # Cylindrical modal coefs for rigid cylinder, but at different distances
        b_N = np.zeros((K, N_order + 1, N_mic), dtype='complex')
        for nm in range(N_mic):
            r = mic_dirs_rad[nm, 1]
            kr = 2 * np.pi * f * r / masp.c

            # Similar to the sph_modal_coefs for the rigid case
            for n in range(N_order+1):
                jn = jv(n, kr)
                jnprime = jvp(n, kR, 1)
                hn = hankel2(n, kr)
                hnprime = h2vp(n, kR, 1)
                b_N[:, n, nm] = np.power(1j, n) * (jn - (jnprime / hnprime) * hn)

    # Avoid NaNs for very high orders, instead of (very) very small values
    b_N[np.isnan(b_N)] = 0.

    # Compute angular-dependent part of the microphone responses
    H_mic = np.zeros((K, N_mic, N_doa), dtype='complex')
    for nd in range(N_doa):
        # Unit vectors of DOAs and microphones
        azi0 = src_dirs_rad[nd]
        azi = mic_dirs_rad[:, 0]
        angle = azi-azi0

        C = np.zeros((N_order + 1, N_mic))
        for n in range(N_order + 1):
            # Jacobi-Anger expansion
            if n == 0:
                C[n, :] = np.ones(angle.shape)
            else:
                C[n, :] = 2 * np.cos(n * angle)
        # Accumulate across orders
        if same_radius:
            H_mic[:, :, nd] = np.matmul(b_N, C)
        else:
            for nm in range(N_mic):
                H_mic[:, nm, nd] = np.matmul(b_N[:, :, nm], C[:, nm])

    # Handle Nyquist for real impulse response
    tempH_mic = H_mic.copy()
    # TODO: in `simulate_sph_array()` it was real, not abs. Why?
    tempH_mic[-1, :] = np.abs(tempH_mic[-1, :])
    # Conjugate ifft and fftshift for causal IR
    h_mic = np.real(np.fft.fftshift(np.fft.ifft(np.append(tempH_mic, np.conj(tempH_mic[-2:0:-1,:]), axis=0), axis=0), axes=0))

    return h_mic, H_mic
