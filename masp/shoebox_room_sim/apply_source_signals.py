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
#   @file   apply_source_signals.py
#   @author Andrés Pérez-López
#   @date   30/07/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import numpy as np
import scipy.signal

from masp.validate_data_types import _validate_ndarray_3D, _validate_ndarray_2D, _validate_ndarray_4D, _validate_list


def apply_source_signals_array(array_rirs, src_sigs):
    """
    Apply room impulse responses from an array of receivers (microphone arrays) to a set of source signals.

    Parameters
    ----------
    array_rirs : List
      RIR for each receiver element. Length = (nRec)
    src_sigs: ndarray
       Matrix containing the source signals. Dimension = (L_sig, nSrc)

    Returns
    -------
    array_sigs : List
        Source signals subjected to the RIRs. Length = (nRec)

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    Each element in `array_rirs` should be a ndarray wit dimension = (L, nMic, nSrc).
    All values of `nSrc` across receivers should be equal, and also match (shape[1]) in `src_sigs`.

    Each element of the algorithm output `array_sigs` contains the rendering of the sources from a different receiver,
    featuring a ndarray with dimension = (L_rir+L_sig-1, nMic).

    TODO: check return values
    """

    _validate_list('array_rirs', array_rirs)
    nRec = len(array_rirs)
    nSrcs = [array_rirs[i].shape[2] for i in range(nRec)]
    assert nSrcs[1:] == nSrcs[:-1]  # Same number of sources for each receiver
    nSrc = nSrcs[0]
    _validate_ndarray_2D('src_sigs', src_sigs, shape1=nSrc)

    array_sigs = [None] * nRec      # Empty list
    src_array_sigs = [None] * nRec  # Empty list
    L_sig = src_sigs.shape[0]

    for nr in range(nRec):
        temprirs = array_rirs[nr]
        L_rir = temprirs.shape[0]
        nMics = temprirs.shape[1]
        tempsigs = np.zeros((L_sig + L_rir - 1, nMics, nSrc))
        for ns in range(nSrc):
            print('Convolving with source signal: Source ' + str(ns) + ' - Receiver ' + str(nr))
            tempsigs[:, :, ns] = scipy.signal.fftconvolve(temprirs[:, :, ns], src_sigs[:, ns, np.newaxis], axes=0)
        src_array_sigs[nr] = tempsigs
        array_sigs[nr] = np.sum(tempsigs, axis=2)

    # return array_sigs, src_array_sigs
    return array_sigs


def apply_source_signals_mic(mic_rirs, src_sigs):
    """
    Apply microphone room impulse responses to a set of source signals.

    Parameters
    ----------
    mic_rirs : ndarray
       Matrix containing the room impulse responses. Dimension = (L_rir, nRec, nSrc)
    src_sigs: ndarray
       Matrix containing the source signals. Dimension = (L_sig, nSrc)

    Returns
    -------
    mic_sigs : ndarray
        Source signals subjected to the RIRs. Dimension = = (L_rir+L_sig-1, nRec)

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    The number of source positions (shape[2]) in `mic_rirs` should match the number of sources (shape[1]) in `src_sigs`.

    """
    L_rir = mic_rirs.shape[0]
    nRec = mic_rirs.shape[1]
    nSrc = mic_rirs.shape[2]
    L_sig = src_sigs.shape[0]
    _validate_ndarray_3D('mic_rirs', mic_rirs)
    _validate_ndarray_2D('src_sigs', src_sigs, shape1=nSrc)

    mic_sigs = np.zeros((L_sig + L_rir - 1, nRec))
    for nr in range(nRec):
        for ns in range(nSrc):
            print('Convolving with source signal: Source ' + str(ns) + ' - Receiver ' + str(nr))
            mic_sigs[:, nr] = mic_sigs[:, nr] + scipy.signal.fftconvolve(mic_rirs[:, nr, ns], src_sigs[:, ns])

    return mic_sigs


def apply_source_signals_sh(sh_rirs, src_sigs):
    """
    Apply spherical harmonic room impulse responses to a set of source signals.

    Parameters
    ----------
    sh_rirs : ndarray
       Matrix containing the room impulse responses. Dimension = (L_rir, nSH, nRec, nSrc)
    src_sigs: ndarray
       Matrix containing the source signals. Dimension = (L_sig, nSrc)

    Returns
    -------
    sh_sigs : ndarray
        Source signals subjected to the RIRs. Dimension = = (L_rir+L_sig-1, nSH, nRec)

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    The number of source positions (shape[3]) in `sh_rirs` should match the number of sources (shape[1]) in `src_sigs`.

    TODO: check return values
    """

    L_rir = sh_rirs.shape[0]
    nSH = sh_rirs.shape[1]
    nRec = sh_rirs.shape[2]
    nSrc = sh_rirs.shape[3]
    L_sig = src_sigs.shape[0]
    _validate_ndarray_4D('mic_rirs', sh_rirs)
    _validate_ndarray_2D('src_sigs', src_sigs, shape1=nSrc)

    sh_sigs = np.zeros((L_sig + L_rir - 1, nSH, nRec))
    src_sh_sigs = np.zeros((L_sig + L_rir - 1, nSH, nRec, nSrc))
    for nr in range(nRec):
        for ns in range(nSrc):
            print('Convolving with source signal: Source ' + str(ns) + ' - Receiver ' + str(nr))
            # Just apply convolution to non-zero channels (skip empty channels due to different sh orders)
            idx_nonzero = np.sum(sh_rirs[:,:, nr, ns],axis=0) != 0
            src_sh_sigs[:, idx_nonzero, nr, ns] = scipy.signal.fftconvolve(sh_rirs[:, idx_nonzero, nr, ns], src_sigs[:, ns, np.newaxis], axes=0)
    sh_sigs = np.sum(src_sh_sigs, axis=3)

    # return sh_sigs, src_sh_sigs
    return sh_sigs

