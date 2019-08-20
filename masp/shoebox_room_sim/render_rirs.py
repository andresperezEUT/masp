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
#   @file   render_rirs.py
#   @author Andrés Pérez-López
#   @date   30/07/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import numpy as np
import scipy.signal
from masp.utils import lagrange

def render_rirs_array(echograms, band_centerfreqs, fs, grids, array_irs):
    """
    TODO
    :param echograms:
    :param band_centerfreqs:
    :param fs:
    :param grids:
    :param array_irs:
    :return:
    """
    raise NotImplementedError



def render_rirs_mic(echograms, band_centerfreqs, fs):
    """

    :param echograms:
    :param band_centerfreqs:
    :param fs:
    :return:
    """

    # echograms: [nSrc, nRec, nBands] dimension
    nSrc = echograms.shape[0]
    nRec = echograms.shape[1]
    nBands = echograms.shape[2]

    # Sample echogram to a specific sampling rate with fractional interpolation
    fractional = False # TODO: ARGUMENT

    # Decide on number of samples for all RIRs
    endtime = 0
    for ns in range(nSrc):
        for nr in range(nRec):
            # TODO: this line gets the last time value of the first echogram.
            #  Is that what we want, or is it the maximum of all last time values?
            # temptime = echograms[ns, nr].time(end)
            temptime = echograms[ns, nr, 0].time[-1]
            if temptime > endtime:
                endtime = temptime

    L_rir = int(np.ceil(endtime * fs))
    if nBands > 1:
        L_fbank = 1000 # TODO: is filterbank lenght hardcoded?
    else:
        L_fbank = 0
    L_tot = L_rir + L_fbank

    # Render responses and apply filterbank to combine different decays at different bands1    rirs = np.zeros((L_tot, maxSH, nRec, nSrc))
    rirs = np.empty((L_tot, nRec, nSrc))
    for ns in range(nSrc):
        for nr in range(nRec):

            print('Rendering echogram: Source ' + str(ns) + ' - Receiver ' + str(nr))

            tempIR = np.zeros((L_rir, nBands))
            for nb in range(nBands):
                tempIR[:, nb] = np.squeeze(render_rirs(echograms[ns, nr, nb], endtime, fs, fractional))
            # TODO: CHECK SQUEEZE...
            print('     Filtering and combining bands')
            rirs[:, nr, ns] = filter_rirs(tempIR, band_centerfreqs, fs)
    return rirs


def render_rirs_sh(echograms, band_centerfreqs, fs):
    """
    TODO
    :param echograms:
    :param band_centerfreqs:
    :param fs:
    :return:
    """

    # echograms: [nSrc, nRec, nBands] dimension
    nSrc = echograms.shape[0]
    nRec = echograms.shape[1]
    nBands = echograms.shape[2]

    # Sample echogram to a specific sampling rate with fractional interpolation
    fractional = False # TODO: ARGUMENT

    # Decide on number of samples for all RIRs
    endtime = 0
    for ns in range(nSrc):
        for nr in range(nRec):
            # TODO: this line gets the last time value of the first echogram.
            #  Is that what we want, or is it the maximum of all last time values?
            # temptime = echograms[ns, nr].time(end)
            temptime = echograms[ns, nr, 0].time[-1]
            if temptime > endtime:
                endtime = temptime
    # # python add
    # endtime = int(np.ceil(endtime))

    L_rir = int(np.ceil(endtime * fs))
    if nBands > 1:
        L_fbank = 1000 # TODO: is filterbank lenght hardcoded?
    else:
        L_fbank = 0
    L_tot = L_rir + L_fbank

    # Find maximum number of SH channels of echograms
    maxSH = 0
    for nr in range(nRec):
        # TODO: should the echograms have different column number as a function of sh?
        #  now that seems broken...
        if np.ndim(echograms[0, nr, 0].value) <= 1:
            tempSH = 1
        else:
            tempSH = np.shape(echograms[0, nr, 0].value)[1] # this should work in the general case (broken now?)
        if tempSH > maxSH:
            maxSH = tempSH

    # Render responses and apply filterbank to combine different decays at different bands1    rirs = np.zeros((L_tot, maxSH, nRec, nSrc))

    # TODO: AGAIN HERE SAME PROBLEM AS ABOVE...
    if np.ndim(echograms[ns, nr, 0].value) <= 1:
        nSH = 1
    else:
        nSH = np.shape(echograms[ns, nr, 0].value)[1]

    rirs = np.empty((L_tot, nSH, nRec, nSrc))
    for ns in range(nSrc):
        for nr in range(nRec):

            print('Rendering echogram: Source ' + str(ns) + ' - Receiver ' + str(nr))

            tempIR = np.zeros((L_rir, nSH, nBands))
            for nb in range(nBands):
                tempIR[:, :, nb] = render_rirs(echograms[ns, nr, nb], endtime, fs, fractional)

            print('     Filtering and combining bands')
            for nh in range(nSH):
                rirs[:, nh, nr, ns] = filter_rirs(np.squeeze(tempIR[:, nh, :]), band_centerfreqs, fs)
    return rirs


def render_rirs(echogram, endtime, fs, fractional):
    """
    TODO
    %RENDER_RIR Samples the echogram to a specified sample rate.
%   echogram:       the echogram structure
%   endtime:        time in secs for desired length of the impulse response
%   fs:             sampling rate
%

    :param echogram:
    :param endtime:
    :param fs:
    :param fractional:
    :return:
    """

    # Number of reflections inside the time limit
    # TODO CHECK
    idx_trans = echogram.time[echogram.time < endtime].size

    if np.ndim(echogram.value) <= 1:
        nChannels = 1
    else:
        nChannels = np.shape(echogram.value)[1]

    if fractional:

        # Get lagrange interpolating filter of order 100 (filter length 101)
        # todo: parametrize order?
        # order = 100
        # h_offset = 50
        # h_idx = np.arange(-(order/2),(order/2)+1).astype(int) # only valid for even orders
        order = 100
        h_offset = 50
        h_idx = np.arange(-(order/2),(order/2)+1).astype(int) # only valid for even orders

        # make a filter table for quick access for quantized fractional samples
        fractions = np.linspace(0,1,101)
        H_frac = lagrange(order, 50 + fractions)

        # Initialise array
        tmp_ir = np.zeros((int(np.ceil(endtime * fs) + (2*h_offset)), nChannels))

        for i in range(idx_trans):
            refl_idx = int(np.floor(echogram.time[i] * fs) + 1)
            refl_frac = np.remainder(echogram.time[i] * fs, 1)
            filter_idx = np.argmin(np.abs(refl_frac - fractions))
            h_frac = H_frac[:, filter_idx]
            tmp_ir[h_offset+refl_idx+h_idx-1,:] += h_frac[:,np.newaxis] * echogram.value[i]

        ir = tmp_ir[h_offset:-h_offset,:]

    else:

        L_ir = int(np.ceil(endtime * fs))
        ir = np.zeros((L_ir, nChannels))
        refl_idx = (np.round(echogram.time[:idx_trans] * fs)).astype(int)
        # filter out exceeding indices
        refl_idx = refl_idx[refl_idx < L_ir]
        ir[refl_idx, :] = echogram.value[:refl_idx.size]

    return ir


def render_quantized(qechogram, endtime, fs, FRACTIONAL):
    '''
    TODO
    %RENDER_HRTF Samples the echogram using HRTFs.
%   qechogram:      the quantized echogram structure array
%   endtime:        time in secs for desired length of the impulse response
%   fs:             sampling rate
%   FRACTIONAL:     flag for using fractional delay filters (1) for non-integer
%                   sample reflection times (more accurate), or just
%                   quantise (0) to the closest sample (faster)
    :param qechogram:
    :param endtime:
    :param fs:
    :param FRACTIONAL:
    :return:
    '''
    raise NotImplementedError


def filter_rirs(rir, f_center, fs):
    """
    TODO
    :param rir:
    :param f_center:
    :param fs:
    :return:
    """

    nBands = rir.shape[1]

    if f_center.size != nBands:
        raise ValueError('The number of bands should match the number of columns in the 2nd dimension of rir')

    if nBands == 1:
        rir_full = rir
    else:
        # Order of filters
        order = 1000
        filters = np.zeros((order + 1, nBands))
        for i in range(nBands):
            if i == 0:
                fl = 30. # TODO: HARDCODED?
                fh = np.sqrt(f_center[i] * f_center[i + 1])
                wl = fl / (fs / 2.)
                wh = fh / (fs / 2.)
                # w = np.array([wl, wh])
                filters[:, i] = scipy.signal.firwin(order+1, [wl, wh], pass_zero='bandpass')
                # filters[:, i] = fir1(order, w, 'bandpass')
            elif i == nBands-1:
                fl = np.sqrt(f_center[i] * f_center[i - 1])
                w = fl / (fs / 2.)
                filters[:, i] = scipy.signal.firwin(order+1, w, pass_zero='highpass')
                # filters[:, i] = fir1(order, w, 'high')
            else:
                fl = np.sqrt(f_center[i] * f_center[i - 1])
                fh = np.sqrt(f_center[i] * f_center[i + 1])
                wl = fl / (fs / 2.)
                wh = fh / (fs / 2.)
                # w = np.array([wl, wh])
                filters[:, i] = scipy.signal.firwin(order + 1, [wl, wh], pass_zero='bandpass')
                # filters[:, i] = fir1(order, w, 'bandpass')

        temp_rir = np.append(rir, np.zeros((order, nBands)), axis=0)
        # rir_filt = fftfilt(filters, temp_rir)
        # rir_filt = scipy.signal.filtfilt(filters, temp_rir)
        rir_filt = scipy.signal.fftconvolve(filters, temp_rir)[:temp_rir.shape[0]]
        rir_full = np.sum(rir_filt, axis=1)
    return rir_full