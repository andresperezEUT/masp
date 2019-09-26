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
import copy

from masp.shoebox_room_sim.echogram import Echogram
from masp.shoebox_room_sim.quantise import get_echo2gridMap, quantise_echogram
from masp.utils import lagrange, C
from masp.validate_data_types import _validate_echogram, _validate_float, _validate_int, _validate_boolean, \
    _validate_ndarray_2D, _validate_ndarray_1D, _validate_echogram_array, _validate_list, \
    _validate_quantised_echogram_array, _validate_ndarray_3D


def render_rirs_array(echograms, band_centerfreqs, fs, grids, array_irs):
    """
    Render the echogram IRs of an array of mic arrays with arbitrary geometries and transfer functions.

    Parameters
    ----------
    echograms : ndarray, dtype = Echogram
        Target echograms. Dimension = (nSrc, nRec, nBands)
    band_centerfreqs : ndarray
        Center frequencies of the filterbank. Dimension = (nBands)
    fs : int
        Target sampling rate
    grids : List
        DoA grid for each receiver. Length = (nRec)
    array_irs : List
        IR of each element of the eceivers. Length = (nRec)

    Returns
    -------
    rirs : List
        RIR for each receiver element. Length = (nRec)

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    For each microphone array (receiver position), we must provide two parameters:
    `grids` contains the angular positions, in azimuth-elevation pairs (in radians),
        of the sampled measured/computed DoAs of the array. ndarray with dimension = (nDoa, C-1).
    `array_irs` is the time-domain IR from each DoA measurement point to each microphone capsule.
        It is therefore a ndarray with dimension = (L1, nMic, nDoa).
    These parameters are independent for each receiver, but `nDoa` must macth within receiver.

    Each of the elements in the algorithm output list `rirs` is a ndarray with dimension = (L2, nMic, nSrc),
    and contains the Room Impulse Response for each capsule/source pair at each receiver (microphone array).

    The highest center frequency must be at most equal to fs/2, in order to avoid aliasing.
    The lowest center frequency must be at least equal to 30 Hz.
    Center frequencies must increase monotonically.

    TODO: expose fractional, L_filterbank as parameter?
    """

    nSrc = echograms.shape[0]
    nRec = echograms.shape[1]
    nBands = echograms.shape[2]

    _validate_echogram_array(echograms)
    _validate_int('fs', fs, positive=True)
    _validate_ndarray_1D('f_center', band_centerfreqs, positive=True, size=nBands, limit=[30,fs/2])
    _validate_list('grids', grids, size=nRec)
    for i in range(nRec):
        _validate_ndarray_2D('grids_'+str(i), grids[i], shape1=C-1)
    _validate_list('array_irs', array_irs, size=nRec)
    for i in range(nRec):
        _validate_ndarray_3D('array_irs_'+str(i), array_irs[i], shape2=grids[i].shape[0])

    # Sample echogram to a specific sampling rate with fractional interpolation
    fractional = True

    # Decide on number of samples for all RIRs
    endtime = 0
    for ns in range(nSrc):
        for nr in range(nRec):
            temptime = echograms[ns, nr, 0].time[-1]
            if temptime > endtime:
                endtime = temptime

    L_rir = int(np.ceil(endtime * fs))
    L_fbank = 1000 if nBands > 1 else 0

    array_rirs = [None] * nRec
    for nr in range(nRec):
        grid_dirs_rad = grids[nr]
        nGrid = np.shape(grid_dirs_rad)[0]
        mic_irs = array_irs[nr]
        L_resp = np.shape(mic_irs)[0]
        nMics = np.shape(mic_irs)[1]
        array_rirs[nr] = np.zeros((L_rir + L_fbank + L_resp - 1, nMics, nSrc))

        for ns in range(nSrc):
            print('Rendering echogram: Source ' + str(ns) + ' - Receiver ' + str(nr) )
            print('      Quantize echograms to receiver grid')
            echo2gridMap = get_echo2gridMap(echograms[ns, nr, 0], grid_dirs_rad)

            tempRIR = np.zeros((L_rir, nGrid, nBands))
            for nb in range(nBands):

                # First step: reflections are quantized to the grid directions
                q_echograms = quantise_echogram(echograms[ns, nr, nb], nGrid, echo2gridMap)
                # Second step: render quantized echograms
                print('      Rendering quantized echograms: Band ' + str(nb))
                tempRIR[:, :, nb], _ = render_quantised(q_echograms, endtime, fs, fractional)

            tempRIR2 = np.zeros((L_rir + L_fbank, nGrid))
            print('      Filtering and combining bands')
            for ng in range(nGrid):
                tempRIR2[:, ng] = filter_rirs(tempRIR[:, ng, :], band_centerfreqs, fs).squeeze()

            # Third step: convolve with directional IRs at grid directions
            idx_nonzero = [i for i in range(tempRIR2.shape[1]) if np.sum(np.power(tempRIR2[:,i], 2)) > 10e-12]   # neglect grid directions with almost no energy
            tempRIR2 = np.row_stack((tempRIR2[:, idx_nonzero], np.zeros((L_resp - 1, len(idx_nonzero)) )))
            for nm in range(nMics):
                tempResp = mic_irs[:, nm, idx_nonzero]
                array_rirs[nr][:, nm, ns] = np.sum(scipy.signal.fftconvolve(tempResp, tempRIR2, axes=0)[:tempRIR2.shape[0], :], axis=1)

    return array_rirs


def render_rirs_mic(echograms, band_centerfreqs, fs):
    """
    Render a mic echogram array into an impulse response matrix.

    Parameters
    ----------
    echograms : ndarray, dtype = Echogram
        Target echograms. Dimension = (nSrc, nRec, nBands)
    band_centerfreqs : ndarray
        Center frequencies of the filterbank. Dimension = (nBands)
    fs : int
        Target sampling rate

    Returns
    -------
    ir : ndarray
        Rendered echograms. Dimension = (M, nRec, nSrc)

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    The highest center frequency must be at most equal to fs/2, in order to avoid aliasing.
    The lowest center frequency must be at least equal to 30 Hz.
    Center frequencies must increase monotonically.

    TODO: expose fractional, L_filterbank as parameter?
    """

    nSrc = echograms.shape[0]
    nRec = echograms.shape[1]
    nBands = echograms.shape[2]
    _validate_echogram_array(echograms)
    _validate_int('fs', fs, positive=True)
    _validate_ndarray_1D('f_center', band_centerfreqs, positive=True, size=nBands, limit=[30,fs/2])

    # Sample echogram to a specific sampling rate with fractional interpolation
    fractional = True

    # Decide on number of samples for all RIRs
    endtime = 0
    for ns in range(nSrc):
        for nr in range(nRec):
            temptime = echograms[ns, nr, 0].time[-1]
            if temptime > endtime:
                endtime = temptime

    L_rir = int(np.ceil(endtime * fs))
    L_fbank = 1000 if nBands > 1 else 0
    L_tot = L_rir + L_fbank

    # Render responses and apply filterbank to combine different decays at different bands
    rirs = np.empty((L_tot, nRec, nSrc))
    for ns in range(nSrc):
        for nr in range(nRec):

            print('Rendering echogram: Source ' + str(ns) + ' - Receiver ' + str(nr))
            tempIR = np.zeros((L_rir, nBands))
            for nb in range(nBands):
                tempIR[:, nb] = np.squeeze(render_rirs(echograms[ns, nr, nb], endtime, fs, fractional))

            print('     Filtering and combining bands')
            rirs[:, nr, ns] = filter_rirs(tempIR, band_centerfreqs, fs).squeeze()

    return rirs


def render_rirs_sh(echograms, band_centerfreqs, fs):
    """
    Render a spherical harmonic echogram array into an impulse response matrix.

    Parameters
    ----------
    echograms : ndarray, dtype = Echogram
        Target echograms. Dimension = (nSrc, nRec, nBands)
    band_centerfreqs : ndarray
        Center frequencies of the filterbank. Dimension = (nBands)
    fs : int
        Target sampling rate

    Returns
    -------
    ir : ndarray
        Rendered echograms. Dimension = (M, maxSH, nRec, nSrc)

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    `maxSH` is the highest spherical harmonic number found in all echograms.
    For any echogram with nSH<maxSH, the channels (nSH...maxSH) will contain only zeros.

    The highest center frequency must be at most equal to fs/2, in order to avoid aliasing.
    The lowest center frequency must be at least equal to 30 Hz.
    Center frequencies must increase monotonically.

    TODO: expose fractional, L_filterbank as parameter?
    """

    # echograms: [nSrc, nRec, nBands] dimension
    nSrc = echograms.shape[0]
    nRec = echograms.shape[1]
    nBands = echograms.shape[2]
    _validate_echogram_array(echograms)
    _validate_int('fs', fs, positive=True)
    _validate_ndarray_1D('f_center', band_centerfreqs, positive=True, size=nBands, limit=[30,fs/2])

    # Sample echogram to a specific sampling rate with fractional interpolation
    fractional = True

    # Decide on number of samples for all RIRs
    endtime = 0
    for ns in range(nSrc):
        for nr in range(nRec):
            temptime = echograms[ns, nr, 0].time[-1]
            if temptime > endtime:
                endtime = temptime

    L_rir = int(np.ceil(endtime * fs))
    L_fbank = 1000 if nBands > 1 else 0
    L_tot = L_rir + L_fbank

    # Find maximum number of SH channels in all echograms
    maxSH = 0
    for nr in range(nRec):
        tempSH = np.shape(echograms[0, nr, 0].value)[1]
        if tempSH > maxSH:
            maxSH = tempSH

    # Render responses and apply filterbank to combine different decays at different bands
    rirs = np.empty((L_tot, maxSH, nRec, nSrc))
    for ns in range(nSrc):
        for nr in range(nRec):

            print('Rendering echogram: Source ' + str(ns) + ' - Receiver ' + str(nr))
            nSH = np.shape(echograms[ns, nr, 0].value)[1]

            tempIR = np.zeros((L_rir, nSH, nBands))
            for nb in range(nBands):
                tempIR[:, :, nb] = render_rirs(echograms[ns, nr, nb], endtime, fs, fractional)

            print('     Filtering and combining bands')
            for nh in range(nSH):
                rirs[:, nh, nr, ns] = filter_rirs(tempIR[:, nh, :], band_centerfreqs, fs).squeeze()
    return rirs




def render_rirs(echogram, endtime, fs, fractional=True):
    """
    Render an echogram into an impulse response.

    Parameters
    ----------
    echogram : Echogram
        Target Echogram.
    endtime : float
        Maximum time of rendered reflections, in seconds.
    fs : int
        Target sampling rate
    fractional : bool, optional
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


def render_quantised(qechogram, endtime, fs, fractional):
    """
    Render a quantised echogram array into a quantised impulse response matrix.

    Parameters
    ----------
    qechograms : ndarray, dtype = QuantisedEchogram
        Target quantised echograms. Dimension = (nDirs).
    endtime : float
        Maximum time of rendered reflections, in seconds.
    fs : int
        Target sampling rate
    fractional : bool, optional
        Use fractional or integer (round) delay. Default to True.

    Returns
    -------
    qIR : ndarray
        Rendered quantised echograms. Dimension = (ceil(endtime * fs), nChannels)
    idx_nonzero : 1D ndarray
        Indices of non-zero elements.

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    TODO: expose fractional as parameter?
    """

    _validate_quantised_echogram_array(qechogram)
    _validate_float('edntime', endtime, positive=True)
    _validate_int('fs', fs, positive=True)
    _validate_boolean('fractional', fractional)

    # Number of grid directions of the quantization
    nDirs = qechogram.size

    # Render echograms
    L_rir = int(np.ceil(endtime * fs))
    qIR = np.zeros((L_rir, nDirs))
    idx_nonzero = []
    for nq in range(nDirs):
        tempgram = Echogram(
            time=qechogram[nq].time,
            value=qechogram[nq].value,
            order=np.zeros((qechogram[nq].time.size, 3), dtype=int),  # whatever to pass the size validation
            coords=np.zeros((qechogram[nq].time.size, 3)))            # whatever to pass the size validation
        # Omit if there are no echoes in the specific one
        if qechogram[nq].isActive:
            idx_nonzero.append(nq)
            # Number of reflections inside the time limit')
            idx_limit = tempgram.time[tempgram.time < endtime].size
            tempgram.time = tempgram.time[:idx_limit+1]
            tempgram.value = tempgram.value[:idx_limit+1]
            tempgram.order = tempgram.order[:idx_limit+1]    # whatever to pass the size validation
            tempgram.coords = tempgram.coords[:idx_limit+1]  # whatever to pass the size validation

            qIR[:, nq] = render_rirs(tempgram, endtime, fs, fractional).squeeze()

    return qIR, np.asarray(idx_nonzero)


def filter_rirs(rir, f_center, fs):
    """
    Apply a filterbank to a given impulse responses.

    Parameters
    ----------
    rir : ndarray
        Impulse responses to be filtered.  Dimension = (L, nBands)
    f_center : ndarray
        Center frequencies of the filterbank. Dimension = (nBands)
    fs : int
        Target sampling rate

    Returns
    -------
    ir : ndarray
        Filtered impulse responses. Dimension = (L+M, 1)

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    Filter operation is implemented with `scipy.signal.firwin`.
    Order of the filters is hardcoded to M = 1000 (length=M+1).

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
