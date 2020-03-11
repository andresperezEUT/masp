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
#   @file   test_script_sh.py
#   @author Andrés Pérez-López
#   @date   24/07/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import numpy as np
from asma import shoebox_room_sim as srs
import time
import librosa
import sys
import scipy.signal
pp = "/Users/andres.perez/source/parametric_spatial_audio_processing"
sys.path.append(pp)
import parametric_spatial_audio_processing as psa
import matplotlib.pyplot as plt
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SETUP

# Room definition
room = np.array([10, 10, 10])

# Desired RT per octave band, and time to truncate the responses
rt60 = np.array([0.5])
nBands = len(rt60)

# Generate octave bands
band_centerfreqs = np.empty(nBands)
band_centerfreqs[0] = 1000

# Absorption for approximately achieving the RT60 above - row per band
abs_wall = srs.find_abs_coeffs_from_rt(room, rt60)[0]

# Critical distance for the room
_, d_critical, _ = srs.room_stats(room, abs_wall)

# Receiver position
rec = np.array([ [5, 5, 5] ])
nRec = rec.shape[0]

# Source positions
src = np.array([ [5, 5, 6] ])
nSrc = src.shape[0]

# SH orders for receivers
rec_orders = np.array([1]) # rec1: first order(4ch), rec2: 3rd order (16ch)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# RUN SIMULATOR

# Echogram
tic = time.time()

maxlim = 1.5 # just stop if the echogram goes beyond that time ( or just set it to max(rt60) )
limits = np.minimum(rt60, maxlim)

# Compute echograms
# abs_echograms, rec_echograms, echograms = srs.compute_echograms_sh(room, src, rec, abs_wall, limits, rec_orders)
abs_echograms = srs.compute_echograms_sh(room, src, rec, abs_wall, limits, rec_orders)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# RENDERING
 
# In this case all the information (e.g. SH directivities) are already
# encoded in the echograms, hence they are rendered directly to discrete RIRs
fs = 48000
sh_rirs = srs.render_rirs_sh(abs_echograms, band_centerfreqs, fs).squeeze()
sh_rirs = sh_rirs * np.sqrt(4*np.pi) * [1, 1./np.sqrt(3), 1./np.sqrt(3), 1./np.sqrt(3)]  # SN3D norm
plt.figure()
plt.plot(sh_rirs)
plt.show()


signal_len_samples = int(np.floor(1. * fs))
signal = np.random.randn(signal_len_samples)

reverberant_signal = np.zeros((signal_len_samples, 4))
for i in range(4):
    reverberant_signal[:,i] = scipy.signal.fftconvolve(signal, sh_rirs[:,i].squeeze())[:signal_len_samples]
x = psa.Signal(reverberant_signal.T, fs, 'acn', 'sn3d')
psa.plot_signal(x,title='waveform')

analysis_window_size = 512
window_overlap = analysis_window_size // 2
fft_size = analysis_window_size
stft = psa.Stft.fromSignal(x,
                           window_size=analysis_window_size,
                           window_overlap=window_overlap,
                           nfft=fft_size)
psa.plot_magnitude_spectrogram(stft,title='magnitude spectrogram')
doa = psa.compute_DOA(stft)
psa.plot_doa(doa,title='doa')

plt.show()

i = psa.compute_intensity_vector(stft)

psa.plot_magnitude_spectrogram(i)
plt.show()