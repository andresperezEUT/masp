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
#   @file   srs_bin_moving.py
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
rec = np.array([3.4, 2.1, 1.7])
nRec = rec.shape[0]

# Surce start
xs0 = np.array([7.4, 2.1, 1.8])
# Source velocity (m/s)
s_vel = 1
# Other point in line
xs1 = np.array([3.4, 3.1, 1.8])
# Duration of source signal (5sec)
tSig = 5
# Block length for RIR interpolation (50msec)
tBlock = 0.05
# number of blocks
nBlocks = int(np.ceil(tSig/tBlock))
# Get source positions for block RIRs
nSrc = nBlocks+1
src = np.zeros((nSrc, 3))
for ns in range(nSrc):
    src[ns,:] = xs0 + (xs1-xs0)/np.sqrt(np.sum(np.power(xs1-xs0, 2)))*ns*tBlock

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# DEFINE DIRECTIONAL IRs FOR RECEIVERS

fs = 48000

# Receiver 3 example: Measured HRTFs
# (MAKE SURE THAT THEY ARE THE SAME SAMPLERATE AS THE REST OR RESAMPLE)
print('Loading measured HRTF responses')
# TODO: WHERE ARE THOSE???????
# load('ownsurround_sim2016.mat','hrtf_dirs','hrtf_mtx')
# grids{1} = hrtf_dirs*pi/180;
# array_irs{1} = hrtf_mtx;
# array_irs{1} = array_irs{1}(250+(1:256),:,:);
# clear hrtf_mtx hrtf_dirs

# TODO: TBC...

