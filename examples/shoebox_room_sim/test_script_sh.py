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
from masp import shoebox_room_sim as srs
import time
import librosa

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SETUP

# Room definition
room = np.array([10.2, 7.1, 3.2])

# Desired RT per octave band, and time to truncate the responses
rt60 = np.array([1., 0.8, 0.7, 0.6, 0.5, 0.4])
nBands = len(rt60)

# Generate octave bands
band_centerfreqs = np.empty(nBands)
band_centerfreqs[0] = 125
for nb in range(1, nBands):
    band_centerfreqs[nb] = 2 * band_centerfreqs[nb-1]

# Absorption for approximately achieving the RT60 above - row per band
abs_wall = srs.find_abs_coeffs_from_rt(room, rt60)[0]

# Critical distance for the room
_, d_critical, _ = srs.room_stats(room, abs_wall)

# Receiver position
rec = np.array([ [4.5, 3.4, 1.5], [2.0, 3.1, 1.4] ])
nRec = rec.shape[0]

# Source positions
src = np.array([ [6.2, 2.0, 1.8], [7.9, 3.3, 1.75], [5.8, 5.0, 1.9] ])
nSrc = src.shape[0]

# SH orders for receivers
rec_orders = np.array([1, 3]) # rec1: first order(4ch), rec2: 3rd order (16ch)

# TODO: not needed?
# % % convert source directions from listener-centric to room-centric
# % [src_coords(:,1), src_coords(:,2), src_coords(:,3)] = sph2cart(src_dirs(:,1)*pi/180, ...
# %     src_dirs(:,2)*pi/180, src_r);
# % src_coords(:,2) = -src_coords(:,2);
# % src = ones(nSrc,1)*rec(1,:) + src_coords;
# % % check sources
# % for n=1:nSrc
# %     if (src(n,1)>room(1))||(src(n,2)>room(2))||(src(n,3)>room(3))
# %         error('Source coordinates out of room boundaries')
# %     end
# % end
#

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
sh_rirs = srs.render_rirs_sh(abs_echograms, band_centerfreqs, fs)

toc = time.time()
print('Elapsed time is ' + str(toc-tic) + 'seconds.')


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# GENERATE SOUND SCENES
# Each source is convolved with the respective mic IR, and summed with
# the rest of the sources to create the microphone mixed signals


sourcepath = '../../data/milk_cow_blues_4src.wav'
src_sigs = librosa.core.load(sourcepath, sr=None, mono=False)[0].T[:,:nSrc]

sh_sigs = srs.apply_source_signals_sh(sh_rirs, src_sigs)