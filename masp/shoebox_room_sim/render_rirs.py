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
from masp.validate_data_types import _validate_echogram, _validate_float, _validate_int, _validate_boolean, \
    _validate_ndarray_2D, _validate_ndarray_1D


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

    TODO: expose fractional as parameter?
    """

    # echograms: [nSrc, nRec, nBands] dimension
    nSrc = echograms.shape[0]
    nRec = echograms.shape[1]
    nBands = echograms.shape[2]

    # Sample echogram to a specific sampling rate with fractional interpolation
    fractional = True

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
            rirs[:, nr, ns] = filter_rirs(tempIR, band_centerfreqs, fs).squeeze()
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




def render_rirs(echogram, endtime, fs, fractional=True):
    """
    Render an echogram into an impulse response.

    Parameters
    ----------
    echogram : Echogram
        Target Echogram.
    endtime: float
        Maximum time of rendered reflections, in seconds.
    fs: int
        Target sampling rate
    fractional: bool, optional
        Use fractional or integer (round) delay. Default to True.

    Returns
    -------
    ir : ndarray
        Rendered echogram. Dimension = (ceil(endtime * fs), nChannels)

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    TODO: expose filter order as parameter?
    """

    _validate_echogram(echogram)
    _validate_float('endtime', endtime, positive=True)
    _validate_int('fs', fs, positive=True)
    _validate_boolean('fractional', fractional)

    nChannels = 1 if np.ndim(echogram.value) <= 1 else np.shape(echogram.value)[1]
    L_ir = int(np.ceil(endtime * fs))
    ir = np.zeros((L_ir, nChannels))
    # Number of reflections inside the time limit
    idx_trans = echogram.time[echogram.time < endtime].size

    if fractional:
        # Get lagrange interpolating filter of order 100 (filter length 101)
        order = 100
        h_offset = 50
        h_idx = np.arange(-(order/2), (order/2)+1).astype(int)  # only valid for even orders

        # Make a filter table for quick access for quantized fractional samples
        fractions = np.linspace(0, 1, 101)
        H_frac = lagrange(order, 50 + fractions)

        # Initialise array
        tmp_ir = np.zeros((int(L_ir + (2*h_offset)), nChannels))

        for i in range(idx_trans):
            refl_idx = int(np.floor(echogram.time[i] * fs) + 1)
            refl_frac = np.remainder(echogram.time[i] * fs, 1)
            filter_idx = np.argmin(np.abs(refl_frac - fractions))
            h_frac = H_frac[:, filter_idx]
            tmp_ir[h_offset+refl_idx+h_idx-1, :] += h_frac[:, np.newaxis] * echogram.value[i]

        ir = tmp_ir[h_offset:-h_offset, :]

    else:
        refl_idx = (np.round(echogram.time[:idx_trans] * fs)).astype(int)
        # Filter out exceeding indices
        refl_idx = refl_idx[refl_idx < L_ir]
        ir[refl_idx, :] = echogram.value[:refl_idx.size]

    return ir


def render_quantized(echogram, endtime, fs, FRACTIONAL):
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
    Apply a filterbank to a given impulse responses.

    Parameters
    ----------
    rir : ndarray
        Impulse responses to be filtered.  Dimension = (L, nBands)
    f_center: ndarray
        Center frequencies of the filterbank. Dimension = (nBands)
    fs: int
        Target sampling rate

    Returns
    -------
    ir : ndarray
        Filtered impulse responses. Dimension = (L+M-1, 1)

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    Filter operation is implemented with `scipy.signal.firwin`.
    Order of the filters is hardcoded to M = 1000.

    The highest center frequency must be at most equal to fs/2, in order to avoid aliasing.
    The lowest center frequency must be at least equal to 30 Hz.
    Center frequencies must increase monotonically.

    TODO: expose filter order, minimum frequency as parameter?
    """

    nBands = rir.shape[1]
    _validate_ndarray_2D('rir', rir)
    _validate_int('fs', fs, positive=True)
    _validate_ndarray_1D('f_center', f_center, positive=True, size=nBands, limit=[30,fs/2])

    if nBands == 1:
        rir_full = rir
    else:
        order = 1000
        filters = np.zeros((order + 1, nBands))
        for i in range(nBands):
            if i == 0:
                fl = 30. 
                fh = np.sqrt(f_center[i] * f_center[i + 1])
                wl = fl / (fs / 2.)
                wh = fh / (fs / 2.)
                filters[:, i] = scipy.signal.firwin(order+1, [wl, wh], pass_zero='bandpass')
            elif i == nBands-1:
                fl = np.sqrt(f_center[i] * f_center[i - 1])
                w = fl / (fs / 2.)
                filters[:, i] = scipy.signal.firwin(order+1, w, pass_zero='highpass')
            else:
                fl = np.sqrt(f_center[i] * f_center[i - 1])
                fh = np.sqrt(f_center[i] * f_center[i + 1])
                wl = fl / (fs / 2.)
                wh = fh / (fs / 2.)
                filters[:, i] = scipy.signal.firwin(order + 1, [wl, wh], pass_zero='bandpass')
        
        temp_rir = np.append(rir, np.zeros((order, nBands)), axis=0)
        rir_filt = scipy.signal.fftconvolve(filters, temp_rir, axes=0)[:temp_rir.shape[0],:]
        rir_full = np.sum(rir_filt, axis=1)[:,np.newaxis]
    
    return rir_full
