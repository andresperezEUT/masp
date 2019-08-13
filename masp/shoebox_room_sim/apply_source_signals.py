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


def apply_source_signals_array(array_rirs, src_sigs):
    """
    TODO
    :param mic_rirs:
    :param src_sigs:
    :return:
    """

    raise NotImplementedError

#     nRec = len(array_rirs)
#     nSrc = size(array_rirs{1}, 3);
#     nSigs = size(src_sigs, 2);
#     if (
#             nSigs < nSrc), error('the number of source signals should be at least as many as the source number in the simulation'); end
#
#     array_sigs = cell(nRec, 1);
#     src_array_sigs = cell(nRec, 1);
#     L_sig = size(src_sigs, 1);
#     for nr=1:nRec
#     temprirs = array_rirs
#     {nr};
#     L_rir = size(temprirs, 1);
#     nMics = size(temprirs, 2);
#     tempsigs = zeros(L_sig + L_rir - 1, nMics, nSrc);
#     for ns=1:nSrc
#
#
# disp(['Convolving with source signal: Source ' num2str(ns) ' - Receiver ' num2str(nr)])
# tempsigs(:,:, ns) = fftfilt(temprirs(:,:, ns), [src_sigs(:, ns); zeros(L_rir - 1, 1)]);
# end
# src_array_sigs
# {nr} = tempsigs;
# array_sigs
# {nr} = sum(tempsigs, 3);
# end


def apply_source_signals_mic(mic_rirs, src_sigs):
    """
    TODO
    :param mic_rirs:
    :param src_sigs:
    :return:
    """
    nRec = mic_rirs.shape[1]
    nSrc = mic_rirs.shape[2]
    L_rir = mic_rirs.shape[0]

    nSigs = src_sigs.shape[1]
    if nSigs < nSrc:
        raise ValueError('The number of source signals should be at least as many as the source number in the simulation')
    L_sig = src_sigs.shape[0]

    mic_sigs = np.zeros((L_sig + L_rir - 1, nRec))
    for nr in range(nRec):
        for ns in range(nSrc):
            print('Convolving with source signal: Source ' + str(ns) + ' - Receiver ' + str(nr))
            mic_sigs[:, nr] = mic_sigs[:, nr] + scipy.signal.fftconvolve(mic_rirs[:, nr, ns], src_sigs[:, ns])

    return mic_sigs


def apply_source_signals_sh(sh_rirs, src_sigs):
    """
    TODO
    :param sh_rirs:
    :param src_sigs:
    :return:
    """

    raise NotImplementedError

    nRec = sh_rirs.shape[2]
    nSrc = sh_rirs.shape[3]
    nSH = sh_rirs.shape[1]
    L_rir = sh_rirs.shape[0]

    nSigs = src_sigs.shape[1]
    if nSigs < nSrc:
        raise ValueError('The number of source signals should be at least as many as the source number in the simulation')
    L_sig = src_sigs.shape[0]

    sh_sigs = np.zeros((L_sig + L_rir - 1, nSH, nRec))
    src_sh_sigs = np.zeros((L_sig + L_rir - 1, nSH, nRec, nSrc))
    for nr in range(nRec):
        for ns in range(nSrc):
            print('Convolving with source signal: Source ' + str(ns) + ' - Receiver ' + str(nr))
                # todo
                # idx_nonzero = find(sum(squeeze(sh_rirs(:,:, nr, ns)))~ = 0)
                # src_sh_sigs(:, idx_nonzero, nr, ns) = fftfilt(sh_rirs(:, idx_nonzero, nr, ns), [src_sigs(:, ns); zeros(L_rir - 1,                                                                                           1)]);

    sh_sigs = np.sum(src_sh_sigs, axis=3)

    return sh_sigs, src_sh_sigs



