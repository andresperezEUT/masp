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
#   @file   test_script_arrays.py
#   @author Andrés Pérez-López
#   @date   30/07/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import numpy as np
from masp import shoebox_room_sim as srs
from masp import array_response_simulator as ars
from masp.utils import get_capsule_positions, c, load_sph_grid, cart2sph, sph2cart
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
rec = np.array([ [4.5, 3.4, 1.5], [2.0, 3.1, 1.4], [3.3, 2.3, 1.6] ])
nRec = rec.shape[0]

# Source positions
src = np.array([ [6.2, 2.0, 1.8], [5.8, 5.0, 1.9] ])
nSrc = rec.shape[0]


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# DEFINE DIRECTIONAL IRs FOR RECEIVERS

fs = 48000

# Receiver 1 example: Eigenmike SMA
print('Simulating EM32 responses for grid')

mic_dirs_rad = get_capsule_positions('eigenmike')
R = mic_dirs_rad[-1,-1]
arrayType = 'rigid'
f_max = 16000
kR_max = 2*np.pi*f_max*R/c
array_order = int(np.ceil(2*kR_max))  # TODO: this formula resembles Daniel's 2006 (Eq. 14), except for the 2 factor. Why?
L_resp = 1024

# Define grid to simulate responses (or that's coming from measurements directly)
grid = load_sph_grid('../../data/N040_M840_Octa.dat')  # 840 points uniformly distributed, cartesian
grid_sph = cart2sph(grid)[:, :-1]
mic_dirs_rad = mic_dirs_rad[:, :-1]
# Compute array IRs
h_eigen = ars.simulate_sph_array(L_resp, mic_dirs_rad, grid_sph, arrayType, R, array_order, fs)


# Receiver 2 example: Uniform circular array of 8 omnis, radius 10cm
print('Simulating 8ch UCA responses for grid')

mic_dirs_deg = np.zeros((8, 3))
mic_dirs_deg[:, 0] = np.arange(0, 360, 360//8)
mic_dirs_deg[: ,1] = np.zeros(np.shape(mic_dirs_deg)[0])
R = 0.1
mic_dirs_deg[: ,2] = R * np.ones(np.shape(mic_dirs_deg)[0])
mic_xyz = sph2cart(mic_dirs_deg)

# Impulse response parameters
L_resp = 256
# Simulate array using get_array_response()
fDirectivity = lambda angle: 1  # Response of omnidirectional microphone
h_uca, _ = ars.get_array_response(grid, mic_xyz, L_resp, fs, mic_dirs=None, fDir_handle=fDirectivity) #  microphone orientation irrelevant in this case


# Receiver 3 example: Measured HRTFs (MAKE SURE THAT THEY ARE THE SAME SAMPLERATE AS THE REST OR RESAMPLE)
print('Loading measured HRTF responses')
# todo: WHERE ARE THOSE?
# load('APolitis_ownsurround_sim2016.mat','hrtf_dirs','hrtf_mtx')
# grids{3} = hrtf_dirs*pi/180;
# array_irs{3} = hrtf_mtx;
# clear hrtf_mtx hrtf_dirs


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# RUN SIMULATOR

# Echogram
tic = time.time()

# Limit the RIR by reflection order or by time-limit
type = 'maxTime'
maxlim = 1.5  # just cut if it's longer than that ( or set to max(rt60) )
limits = np.minimum(rt60, maxlim)

# Compute echograms
abs_echograms = srs.compute_echograms_array(room, src, rec, abs_wall, limits)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# RENDERING

# In the array case, for each receiver position, a grid of directions
# should be provided, along with the corresponding multichannel IRs for the
# elements of the array at the grid directions.

# TODO: use different data, once we got it...
grids = [grid_sph, grid_sph, grid_sph]
array_irs = [h_uca, h_uca, h_uca]

array_rirs = srs.render_rirs_array(abs_echograms, band_centerfreqs, fs, grids, array_irs)

toc = time.time()
print('Elapsed time is ' + str(toc-tic) + 'seconds.')


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# GENERATE SOUND SCENES

sourcepath = '../../data/milk_cow_blues_4src.wav'
src_sigs = librosa.core.load(sourcepath, sr=None, mono=False)[:3].T[:,:nSrc]

mic_sigs = srs.apply_source_signals_array(array_rirs, src_sigs)

