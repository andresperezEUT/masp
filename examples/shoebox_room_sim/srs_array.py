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
abs_wall = srs.findAbsCoeffsFromRT(room, rt60)[0]

# Critical distance for the room
_, d_critical, _ = srs.room_stats(room, abs_wall)

# Receiver position
rec = np.array([ [4.5, 3.4, 1.5], [2.0, 3.1, 1.4], [3.3, 2.3, 1.6] ])
nRec = rec.shape[0]

# Source positions
src = np.array([ [6.2, 2.0, 1.8], [5.8, 5.0, 1.9] ])
nSrc = rec.shape[0]

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


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# DEFINE DIRECTIONAL IRs FOR RECEIVERS

fs = 48000

# Receiver 1 example: Eigenmike SMA
print('Simulating EM32 responses for grid')
# TODO: DECLARE THIS DATA SOMEWHERE ELSE TO BE REUSED
mic_dirs_deg = np.array([   [0, 32, 0, 328, 0, 45, 69, 45, 0, 315, 291, 315, 91, 90, 90, 89, 180, 212, 180, 148, 180, 225, 249, 225, 180, 135, 111, 135, 269, 270, 270, 271],
                            [21, 0, -21, 0, 58, 35, 0, -35, -58, -35, 0, 35, 69, 32, -31, -69, 21, 0, -21, 0, 58, 35, 0, -35, -58, -35, 0, 35, 69, 32, -32, -69] ])
mic_dirs_rad = mic_dirs_deg*np.pi/180
R = 0.042 #m
arrayType = 'rigid'
c = 343
f_max = 16000
kR_max = 2*np.pi*f_max*R/c
array_order = np.ceil(2*kR_max)
L_resp = 1024

# Define grid to simulate responses (or that's coming from measurements directly)
# TODO: SPHERICAL LIBRARY OR SOMETHING...
# grid = loadSphGrid('N040_M840_Octa.dat'); # 840 points uniformly distributed
# grids{1} = grid.aziElev;

# todo: continue...