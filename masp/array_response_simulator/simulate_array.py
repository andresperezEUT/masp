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
#   @file   simulate_array.py
#   @author Andrés Pérez-López
#   @date   22/08/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import numpy as np


def simulate_cyl_array():
    raise NotImplementedError

def simulate_sph_array(N_filt, mic_dirs_rad, src_dirs_rad, arrayType, R, N_order, fs, dirCoeff=None):
    """
    %SIMULATESPHARRAY Simulate the impulse responses of a spherical array.
    %
    %   SIMULATESPHARRAY computes the impulse responses of the microphones of a
    %   spherical microphone array for the given directions of incident plane
    %   waves. The array type can be either 'open' for omnidirectional
    %   microphones in an open setup, or 'rigid' for omnidirectional
    %   microphones mounted on a sphere.
    %
    %   mic_dirs_rad :  [mic_azi1 mic_elev1; ...] directions of microphone
    %                   capsules
    %
    %   src_dirs_rad :  [src_azi1 src_elev1; ...] DOAs of plane waves for
    %                   evaluation
    %
    %   arrayType    :  'open', 'rigid', or 'directional'
    %   R            :  array radius in m
    %   N_order      :  maximum order of spherical approximation
    %   fs           :  sample rate
    %   dirCoeff     :  if array consists of directional microphones, then
    %                   dirCoeff is the directivity coefficient from 0 to 1 (1
    %                   for omni, 0.5 cardioid, 0 dipole)
    %
    %
    % H_mic:    Output the frequency responses for M directions for the N
    %           microphones. Only half the FFT spectrum is returned (up to Nyquist).
    %           Last dimension is the DOA.
    % h_mic:    Output the IRs for M directions for the N microphones. Last
    %           dimension is the DOA.
    %

    TODO
    """

    if dirCoeff is None:
        dirCoeff = []

    # Compute the frequency-dependent part of the microphone responses(radial dependence)

    f = np.arange(N_filt/2) * fs / N_filt
    c = 343.
    kR = 2*np.pi*f*R/c
    b_N = sph_modal_coefs(N_order, kR, arrayType, dirCoeff)
