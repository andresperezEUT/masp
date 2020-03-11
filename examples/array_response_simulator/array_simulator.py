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
#   @file   array_simulator.py
#   @author Andrés Pérez-López
#   @date   30/08/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# This is a collection of python methods for simulation of array responses of
#
#   a) directional sensors and,
#
#   b) sensors mounted on, or at a distance from, a rigid spherical/cylindrical
#   scatterer.
#
# The computation of their frequency and impulse responses is
# based on the theoretical expansion of a scalar incident plane wave field to
# a series of wavenumber-dependent Bessel-family functions and
# direction-dependent Fourier or Legendre functions.
#
# Alternatively, a function for arbitrary open arrays of directional microphones
# is included, not based on the expansion but directly on the steering vector
# formula of inter-sensor delays and sensor gains for user-defined directional
# patterns (e.g. arrays of cardioid microphones).
#
# For more information on the expansions, you can have a look on [ref.1]
# and [ref.2]
#
# and for example on
# <http://en.wikipedia.org/wiki/Plane_wave_expansion> and
# <http://en.wikipedia.org/wiki/Jacobi-Anger_expansion>
#
# For any questions, comments, corrections, or general feedback, please
# contact archontis.politis@aalto.fi
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# REFERENCES
#
# [1] Earl G. Williams, "Fourier Acoustics: Sound Radiation and Nearfield
#     Acoustical Holography", Academic Press, 1999
#
# [2] Heinz Teutsch, "Modal Array Signal Processing: Principles and
#     Applications of Acoustic Wavefield Decomposition", Springer, 2007
#
# [3] Elko, W. Gary,  "Differential Microphone Arrays",
#     In Y. Huang & J. Benesty (Eds.), Audio signal processing for
#     next-generation multimedia communication systems (pp. 11?65).
#     Springer, 2004
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import numpy as np
from asma.utils import sph2cart, get_capsule_positions
from asma import array_response_simulator as asr
import matplotlib.pyplot as plt
plt.set_cmap('jet')


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  ARRAYS OF OMNIDIRECTIONAL MICROPHONES
#
# This is the simplest case and the array responses can be computed
# either by `simulate_cyl_array()` for 2D, `simulate_sph_array()` for 3D, or
# `get_array_response()` for both. Since `get_array_response()` does not involve a
# spherical or cylindrical expansion, it is the easiest to use in this case.

# A uniform circular array of 5 elements at 10cm radius
mic_azi = np.arange(0,360,360/5)*np.pi/180
mic_elev = np.zeros(np.size(mic_azi))
R = 0.1
R_mic = sph2cart(np.column_stack((mic_azi, mic_elev, R * np.ones(np.size(mic_azi)))))

# Plane wave direcitons to evaluate response
doa_azi = np.arange(0,360,5)*np.pi/180
doa_elev = np.zeros(np.size(doa_azi))
r = np.ones(np.size(doa_azi))
U_doa = sph2cart(np.column_stack((doa_azi, doa_elev, r)))

# Impulse response parameters
fs = 48000
Lfilt = 200

# Simulate array using `get_array_response()`
fDirectivity = lambda angle: 1 # response of omnidirectional microphone
h_mic1, H_mic1 = asr.get_array_response(U_doa, R_mic,
                                        Lfilt, fs,
                                        mic_dirs=None,  # microphone orientation irrelevant in this case
                                        fDir_handle=fDirectivity)

# Simulate array using 2D propagation and cylindrical expansion
N_order = 40  # order of expansion
h_mic2, H_mic2 = asr.simulate_cyl_array(Lfilt, mic_azi, doa_azi, 'open', R, N_order, fs)

# Simulate array using 3D propagation and spherical expansion
h_mic3, H_mic3 = asr.simulate_sph_array(Lfilt,
                                        np.column_stack((mic_azi, mic_elev)),
                                        np.column_stack((doa_azi, doa_elev)),
                                        'open', R, N_order, fs)

# Plots
Nfft = Lfilt
f = np.arange(Nfft//2+1) * fs / Nfft

plt.figure(figsize=plt.figaspect(1/2.))
plt.subplot(2,3,1)
plt.pcolormesh(doa_azi*180/np.pi, f, 20*np.log10(np.abs(H_mic1[:,0,:])), vmin=-20, vmax=10)
plt.colorbar()
plt.ylabel('Frequency (Hz)')
plt.semilogy()
plt.ylim([50, 20000])
plt.xlim([0, 355])
plt.xlabel('Azimuth (radians)')
plt.title('Direct evaluation, magnitude')

plt.subplot(2,3,4)
plt.pcolormesh(doa_azi*180/np.pi, f, np.angle(H_mic1[:,0,:]))
plt.colorbar()
plt.ylabel('Frequency (Hz)')
plt.semilogy()
plt.ylim([50, 20000])
plt.xlim([0, 355])
plt.xlabel('Azimuth (radians)')
plt.title('Direct evluation, phase')

plt.subplot(2,3,2)
plt.pcolormesh(doa_azi*180/np.pi, f, 20*np.log10(np.abs(H_mic2[:,0,:])), vmin=-20, vmax=10)
plt.colorbar()
plt.ylabel('Frequency (Hz)')
plt.semilogy()
plt.ylim([50, 20000])
plt.xlim([0, 355])
plt.xlabel('Azimuth (radians)')
plt.title('Direct evaluation, magnitude')

plt.subplot(2,3,5)
plt.pcolormesh(doa_azi*180/np.pi, f, np.angle(H_mic2[:,0,:]))
plt.colorbar()
plt.ylabel('Frequency (Hz)')
plt.semilogy()
plt.ylim([50, 20000])
plt.xlim([0, 355])
plt.xlabel('Azimuth (radians)')
plt.title('Direct evluation, phase')

plt.subplot(2,3,3)
plt.pcolormesh(doa_azi*180/np.pi, f, 20*np.log10(np.abs(H_mic3[:,0,:])), vmin=-20, vmax=10)
plt.colorbar()
plt.ylabel('Frequency (Hz)')
plt.semilogy()
plt.ylim([50, 20000])
plt.xlim([0, 355])
plt.xlabel('Azimuth (radians)')
plt.title('Direct evaluation, magnitude')

plt.subplot(2,3,6)
plt.pcolormesh(doa_azi*180/np.pi, f, np.angle(H_mic3[:,0,:]))
plt.colorbar()
plt.ylabel('Frequency (Hz)')
plt.semilogy()
plt.ylim([50, 20000])
plt.xlim([0, 355])
plt.xlabel('Azimuth (radians)')
plt.title('Direct evluation, phase')
plt.show()

# Plot the array impulse responses for all microphones and doa at 45deg
plt.figure()
plt.plot(np.squeeze(h_mic1[:,:,45//5]))
plt.suptitle('IRs, 5 element UCA of omnidirectional sensors, DOA at 45deg')
plt.show()

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ARRAY OF DIRECTIONAL MICROPHONES
#
# Simulate the same case but for first-order directional microphones
# instead of omnidirectional ones. These can be evaluated using
# `get_array_response()` or by spherical expansion using `simulate_sph_array()`.
# At the moment circular expansion of directional microphones is not
# included, and since there are two ways already to compute, it's probably
# not needed.

# First-order directivity coefficient: a + (1-a)cos(theta).
# For a = 1, same as the previous omnidirectional case above
a = 0.5  # cardioid response

# Simulate array using `get_array_response()`

fDirectivity = lambda angle: a + (1 - a) * np.cos(angle) # response for cardioid microphone
h_mic1, H_mic1 = asr.get_array_response(U_doa, R_mic,
                                        Lfilt, fs,
                                        mic_dirs=None,  # microphone orientation radial
                                        fDir_handle=fDirectivity)


# simulate array using spherical expansion
arrayType = 'directional'
h_mic2, H_mic2 = asr.simulate_sph_array(Lfilt,
                                        np.column_stack((mic_azi, mic_elev)),
                                        np.column_stack((doa_azi, doa_elev)),
                                        arrayType, R, N_order, fs, dirCoef=a)

# Plots
plt.figure(figsize=plt.figaspect(1/2.))
plt.subplot(2,2,1)
plt.pcolormesh(doa_azi*180/np.pi, f, 20*np.log10(np.abs(H_mic1[:,0,:])), vmin=-20, vmax=10)
plt.colorbar()
plt.ylabel('Frequency (Hz)')
plt.semilogy()
plt.ylim([50, 20000])
plt.xlim([0, 355])
plt.xlabel('Azimuth (radians)')
plt.title('Direct evaluation, magnitude')

plt.subplot(2,2,3)
plt.pcolormesh(doa_azi*180/np.pi, f, np.angle(H_mic1[:,0,:]))
plt.colorbar()
plt.ylabel('Frequency (Hz)')
plt.semilogy()
plt.ylim([50, 20000])
plt.xlim([0, 355])
plt.xlabel('Azimuth (radians)')
plt.title('Direct evluation, phase')

plt.subplot(2,2,2)
plt.pcolormesh(doa_azi*180/np.pi, f, 20*np.log10(np.abs(H_mic2[:,0,:])), vmin=-20, vmax=10)
plt.colorbar()
plt.ylabel('Frequency (Hz)')
plt.semilogy()
plt.ylim([50, 20000])
plt.xlim([0, 355])
plt.xlabel('Azimuth (radians)')
plt.title('Spherical expansion, magnitude')

plt.subplot(2,2,4)
plt.pcolormesh(doa_azi*180/np.pi, f, np.angle(H_mic2[:,0,:]))
plt.colorbar()
plt.ylabel('Frequency (Hz)')
plt.semilogy()
plt.ylim([50, 20000])
plt.xlim([0, 355])
plt.xlabel('Azimuth (radians)')
plt.title('Spherical expansion, phase')

# Plot the array impulse responses for all microphones and doa at 45deg
plt.figure()
plt.plot(np.squeeze(h_mic1[:,:,45//5]))
plt.suptitle('IRs, 5 element UCA of cardioid sensors, DOA at 45deg')
plt.show()

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ARRAY OF HIGHER ORDER DIRECTIONAL MICROPHONES
#
# Higher-order directional patterns can be simulated using `get_array_response()`.
# Appropriate functions should be defined, and different patterns can be
# passed for each element of the array. Not that this approach will
# ignore the frequency dependency that occurs in practice with the
# acoustical processing generating these higher-order patterns. For
# proper response modeling in such cases, one can simulate the response
# of the elements of the subarray, and apply the associated processing for
# pattern generation to the resulting IRs (see e.g. [ref.3]).

# Define an array of two differential patterns, one 2nd-order supercardioid
# and a 3rd-order hypercardioid, oriented at front and at 10cm distance
# between (weights from [ref.3])
a2_scard = np.asarray([(3 - np.sqrt(7))/4, (3*np.sqrt(7) - 7)/2, (15 - (5*np.sqrt(7)))/4])
a3_hcard = np.asarray([-3/32, -15/32, 15/32, 35/32])

fDir1 = lambda angle: a2_scard[0] + a2_scard[1]*np.cos(angle) + a2_scard[2]*np.power(np.cos(angle), 2)
fDir2 = lambda angle: a3_hcard[0] + a3_hcard[1]*np.cos(angle) + \
                      a3_hcard[2]*np.power(np.cos(angle), 2) + a3_hcard[3]*np.power(np.cos(angle), 3)

fdirHandles = np.asarray([fDir1, fDir2])
R_mic = 0.5 * np.asarray([ [0, 1, 0], [0, -1, 0] ])
U_orient = np.asarray([ [1, 0, 0], [1, 0, 0] ])

# Simulate array using getArrayResponse()
_, H_mic = asr.get_array_response(U_doa, R_mic, Lfilt, fs, U_orient, fdirHandles)

# Plots
plt.figure(figsize=plt.figaspect(1/2.))
plt.subplot(2,2,1)
plt.pcolormesh(doa_azi*180/np.pi, f, 20*np.log10(np.abs(H_mic[:,0,:])), vmin=-40, vmax=10)
plt.colorbar()
plt.ylabel('Frequency (Hz)')
plt.semilogy()
plt.ylim([50, 20000])
plt.xlim([0, 355])
plt.xlabel('Azimuth (radians)')
plt.title('2nd-order supercardioid, magnitude')

plt.subplot(2,2,3)
plt.pcolormesh(doa_azi*180/np.pi, f, np.angle(H_mic[:,0,:]))
plt.colorbar()
plt.ylabel('Frequency (Hz)')
plt.semilogy()
plt.ylim([50, 20000])
plt.xlim([0, 355])
plt.xlabel('Azimuth (radians)')
plt.title('2nd-order supercardioid, phase')

plt.subplot(2,2,2)
plt.pcolormesh(doa_azi*180/np.pi, f, 20*np.log10(np.abs(H_mic[:,1,:])), vmin=-40, vmax=10)
plt.colorbar()
plt.ylabel('Frequency (Hz)')
plt.semilogy()
plt.ylim([50, 20000])
plt.xlim([0, 355])
plt.xlabel('Azimuth (radians)')
plt.title('3rd-order hypercardioid, magnitude')

plt.subplot(2,2,4)
plt.pcolormesh(doa_azi*180/np.pi, f, np.angle(H_mic[:,1,:]))
plt.colorbar()
plt.ylabel('Frequency (Hz)')
plt.semilogy()
plt.ylim([50, 20000])
plt.xlim([0, 355])
plt.xlabel('Azimuth (radians)')
plt.title('3rd-order hypercardioid, phase')
plt.show()

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ARRAY OF MICROPHONES ON A CYLINDRICAL/SPHERICAL SCATTERER
#
# Microphones mounted on a spherical scatterer (for 3D processing), or on
# a cylindrical scatterer (for 2D processing) exhibit some advantages over
# open arrays of omnidirectional microphones, mainly due to the inherent
# directionality of the scattering body, and the fact that they are more
# suitable for working in the spatial Fourier transform domain (also
# known as eigenbeam processing, or phase-mode processing in the
# literature) in which case the open array exhibits resonant frequencies
# for which the transformed signals vanish.
#
# The response of a plane wave scattered by these fundamental geometries
# is known analytically, and it involves the special functions and their
# derivatives, included in the library.

# Simulate a single microphone mounted on a cylinder and on a sphere of radius 5cm
R = 0.05
arrayType = 'rigid'
N_order = 40  # order of expansion
mic_azi = np.asarray([0])
mic_dir = np.asarray([[0, 0]])

# Compute simulated impulse responses
doa_azi = np.arange(0,360,5)*np.pi/180
doa_dirs = np.column_stack((doa_azi, np.zeros(np.size(doa_azi))))
_, H_mic_cyl = asr.simulate_cyl_array(Lfilt, mic_azi, doa_azi, arrayType, R, N_order, fs)
_, H_mic_sph = asr.simulate_sph_array(Lfilt, mic_dir, doa_dirs, arrayType, R, N_order, fs)

# Plots
plt.figure(figsize=plt.figaspect(1/2.))
plt.subplot(1,2,1)
plt.pcolormesh(doa_azi*180/np.pi, f, 20*np.log10(np.abs(H_mic_sph[:,0,:])), vmin=-20, vmax=10)
plt.colorbar()
plt.ylabel('Frequency (Hz)')
plt.semilogy()
plt.ylim([50, 20000])
plt.xlim([0, 355])
plt.xlabel('Azimuth (radians)')
plt.title('Spherical Array, mic_dir = [0 0], R = 5cm')

plt.subplot(1,2,2)
plt.pcolormesh(doa_azi*180/np.pi, f, 20*np.log10(np.abs(H_mic_cyl[:,0,:])), vmin=-20, vmax=10)
plt.colorbar()
plt.ylabel('Frequency (Hz)')
plt.semilogy()
plt.ylim([50, 20000])
plt.xlim([0, 355])
plt.xlabel('Azimuth (radians)')
plt.title('Cylindrical Array, mic_dir = [0 0], R = 5cm')
plt.show()

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ARRAY OF MICROPHONES AT A DISTANCE FROM A CYLINDRICAL/SPHERICAL SCATTERER
#
# Apart from the functions that evaluate the responses on the surface of a
# scatterer, functions are included that evaluate the response at some
# distance away from the scatterer, for which the effect gradually
# diminishes with distance. This case may be useful for example in simulating
# arrays that place microphones at different radii for broadband eigenbeam
# processing, or for simulating a hearing aid array from some small
# distance of a spherical head.

# Simulate a microphone mounted on a cylinder and on a sphere of radius
# 5cm, and a second microphone at 5cm distance from the scatterer
R = 0.05
N_order = 40  # order of expansion
mic_azi = 0
mic_elev = 0
mic_dirs_sph = np.asarray([ [mic_azi, mic_elev, R], [mic_azi, mic_elev, 2*R] ])
mic_dirs_cyl = np.asarray([ [mic_azi, R], [mic_azi, 2*R] ])

# Compute simulated impulse responses
doa_azi = np.arange(0,360,5)*np.pi/180
doa_dirs = np.column_stack((doa_azi, np.zeros(np.size(doa_azi))))

_, H_mic_cyl = asr.cylindrical_scatterer(mic_dirs_cyl, doa_azi, R, N_order, Lfilt, fs)
_, H_mic_sph = asr.spherical_scatterer(mic_dirs_sph, doa_dirs, R, N_order, Lfilt, fs)


plt.figure(figsize=plt.figaspect(1/2.))
plt.subplot(2,2,1)
plt.pcolormesh(doa_azi*180/np.pi, f, 20*np.log10(np.abs(H_mic_sph[:,0,:])), vmin=-20, vmax=10)
plt.colorbar()
plt.ylabel('Frequency (Hz)')
plt.semilogy()
plt.ylim([50, 20000])
plt.xlim([0, 355])
plt.xlabel('Azimuth (radians)')
plt.title('Spherical Array, mic_dir = [0 0], R_scat = 5cm, r_mic = 5cm')

plt.subplot(2,2,2)
plt.pcolormesh(doa_azi*180/np.pi, f, 20*np.log10(np.abs(H_mic_cyl[:,0,:])), vmin=-20, vmax=10)
plt.colorbar()
plt.ylabel('Frequency (Hz)')
plt.semilogy()
plt.ylim([50, 20000])
plt.xlim([0, 355])
plt.xlabel('Azimuth (radians)')
plt.title('Cylindrical Array, mic_dir = [0 0], R_scat = 5cm, r_mic = 5cm')

plt.subplot(2,2,3)
plt.pcolormesh(doa_azi*180/np.pi, f, 20*np.log10(np.abs(H_mic_sph[:,1,:])), vmin=-20, vmax=10)
plt.colorbar()
plt.ylabel('Frequency (Hz)')
plt.semilogy()
plt.ylim([50, 20000])
plt.xlim([0, 355])
plt.xlabel('Azimuth (radians)')
plt.title('Spherical Array, mic_dir = [0 0], R_scat = 5cm, r_mic = 10cm')

plt.subplot(2,2,4)
plt.pcolormesh(doa_azi*180/np.pi, f, 20*np.log10(np.abs(H_mic_cyl[:,1,:])), vmin=-20, vmax=10)
plt.colorbar()
plt.ylabel('Frequency (Hz)')
plt.semilogy()
plt.ylim([50, 20000])
plt.xlim([0, 355])
plt.xlabel('Azimuth (radians)')
plt.title('Spherical Array, mic_dir = [0 0], R_scat = 5cm, r_mic = 10cm')
plt.show()


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3D UNIFORM SPHERICAL ARRAY EXAMPLE
#
# Simulate a uniform spherical array with the specifications of the Eigenmike
# array, suitable for up to 4th-order eigenbeam processing

# Eigenmike positions
mic_dirs_rad = get_capsule_positions('eigenmike')

# Type and order of approximation
arrayType = 'rigid'
N_order = 30

# Obtain responses for front, side and up DOAs
src_dirs_rad = np.asarray([ [0, 0], [np.pi, 0], [0, np.pi/2] ])
[h_mic, H_mic] = asr.simulate_sph_array(Lfilt, mic_dirs_rad[:, :-1], src_dirs_rad, arrayType, R, N_order, fs);


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# REFERENCES
#
# [1] Earl G. Williams, "Fourier Acoustics: Sound Radiation and Nearfield
#     Acoustical Holography", Academic Press, 1999
#
# [2] Heinz Teutsch, "Modal Array Signal Processing: Principles and
#     Applications of Acoustic Wavefield Decomposition", Springer, 2007
#
# [3] Elko, W. Gary,  "Differential Microphone Arrays",
#     In Y. Huang & J. Benesty (Eds.), Audio signal processing for
#     next-generation multimedia communication systems (pp. 11?65).
#     Springer, 2004
