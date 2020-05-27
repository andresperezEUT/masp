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
#   @file   acoustics.py
#   @author Andrés Pérez-López
#   @date   30/07/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import numpy as np
import scipy.optimize

# from masp.shoebox_room_sim.validate_data_types import _validate_room
# from masp.shoebox_room_sim.validate_data_types import _validate_rt60_target
# from masp.shoebox_room_sim.validate_data_types import _validate_abs_wall_ratios
# from masp.shoebox_room_sim.validate_data_types import _validate_alpha
# from masp.shoebox_room_sim.validate_data_types import _validate_alpha_walls_per_band

from masp.validate_data_types import _validate_ndarray_1D, _validate_number, _validate_ndarray_2D
from masp.utils import C, c

def find_abs_coeffs_from_rt(room, rt60_target, abs_wall_ratios=None):
    """
    Compute wall absorption coefficients per frequency band and wall.

    Parameters
    ----------
    room : ndarray
        Room dimensions in cartesian coordinates. Dimension = (3) [x, y, z].
    rt60_target : ndarray
        Target reverberation time. Dimension = (nBands).
    abs_wall_ratios : ndarray, optional
        Wall absorption coefficient ratios. Dimension = (6).

    Returns
    -------
    alpha_walls : ndarray
        Wall absorption coefficients . Dimension = (nBands, 6).
    rt60_true : ndarray
        RT60 time computed from result. Dimension = (nBands).

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    nBands will be determined by the length of rt60_target.

    If `abs_wall_ratios` is not specified, no wall absorption is applied.

    abs_wall_ratios are expected to be normalized to 1.
    The method will automatically normalize them, in case.

    """

    # Validate arguments
    _validate_ndarray_1D('room', room, size=3, positive=True)
    _validate_ndarray_1D('rt60_target', rt60_target, positive=True)
    if abs_wall_ratios is not None:
        _validate_ndarray_1D('abs_wall_ratios', abs_wall_ratios, size=6, positive=True)


    # Default wall absorption
    if abs_wall_ratios is None:
        abs_wall_ratios = np.ones(6)

    # Normalize
    abs_wall_ratios = abs_wall_ratios / np.max(abs_wall_ratios)

    nBands = len(rt60_target)
    rt60_true = np.zeros(nBands)
    alpha_walls = np.zeros((nBands,6))

    for nb in range(nBands):
        rt60 = rt60_target[nb]
        fmin = lambda alpha: np.abs(rt60 - get_rt_sabine(alpha, room, abs_wall_ratios))
        alpha = scipy.optimize.fmin(func=fmin, x0=0.0001, disp=False)
        rt60_true[nb] = rt60 + fmin(alpha)
        alpha_walls[nb,:] = alpha * abs_wall_ratios

    return alpha_walls, rt60_true



def get_rt_sabine(alpha, room, abs_wall_ratios):
    """
    Estimate RT60 through Sabine's method.

    Parameters
    ----------
    alpha: int, float or 1-D ndarray
        Absorption coefficient.
    room : ndarray
        Room dimensions in cartesian coordinates. Dimension = (3) [x, y, z].
    abs_wall_ratios : ndarray
        Wall absorption coefficients, in the range [0,1]. Dimension = (6).

    Returns
    -------
    rt60 : float
        Estimated reverberation time.

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    As opposed to `find_abs_coeffs_from_rt()`, `abs_wall_ratios` must be explicit.

    `abs_wall_ratios` must have all values in the range [0,1].

    """

    # Validate arguments
    _validate_number('alpha', alpha)
    _validate_ndarray_1D('room', room, size=C, positive=True)
    _validate_ndarray_1D('abs_wall_ratios', abs_wall_ratios, size=2*C, norm=True)

    l, w, h = room
    V = l*w*h   # room volume
    Stot = 2 * ( (l * w) + (l * h) + (w * h) ) # room area

    alpha_walls = alpha * abs_wall_ratios
    a_x = alpha_walls[[0,1]]
    a_y = alpha_walls[[2,3]]
    a_z = alpha_walls[[4,5]]

    # Mean absorption
    a_mean = np.sum( (w * h * a_x) + (l * h * a_y) + (l * w * a_z) ) / Stot
    rt60 = (55.25 * V) / ( c * Stot * a_mean )

    return rt60


def room_stats(room, abs_wall, verbose=True):
    """
    Estimate RT60 through Sabine's method.

    Parameters
    ----------
    room : ndarray
        Room dimensions in cartesian coordinates. Dimension = (3) [x, y, z].
    abs_wall : ndarray
        Wall absorption coefficients per band. Dimension = (nBands, 6)
    verbose: bool, optional
        Display room stats. Default to False.

    Returns
    -------
    rt60 : float
        Estimated reverberation time.
    d_critical: float
        Estimated critical distance.
    d_mfpath: float
        Estimated mean free path.

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    As opposed to `find_abs_coeffs_from_rt()`, `abs_wall_ratios` must be explicit.

    `alpha` and `abs_wall_ratios` must have all values in the range [0,1].

    """

    _validate_ndarray_1D('room', room, size=C, positive=True)
    _validate_ndarray_2D('abs_wall', abs_wall, shape1=2*C, norm=True)

    l, w, h = room
    V = l*w*h   # room volume
    Stot = 2 * ( (l * w) + (l * h) + (w * h) ) # room area

    # Analyse in frequency bands
    nBands = abs_wall.shape[0]
    a_mean = np.empty(nBands)
    rt60_sabine = np.empty(nBands)
    for m in range(nBands):
        a_x = abs_wall[m,0:2]
        a_y = abs_wall[m,2:4]
        a_z = abs_wall[m,4:6]
        # mean absorption
        a_mean[m] = sum( (w * h * a_x) + (l * h * a_y) + (l * w * a_z)) / Stot
        rt60_sabine[m] = (55.25 * V) / (c * Stot * a_mean[m])

    d_critical = 0.1 * np.sqrt(V / (np.pi * rt60_sabine))
    d_mfpath = 4 * V / Stot

    if verbose:
        print('Room dimensions (m)          ' + str(l) + 'x' + str(w) + 'x' + str(h))
        print('Room volume (m^3)            ' + str(V))
        print('Mean absorption coeff        ' + str(a_mean))
        print('Sabine Rev. Time 60dB (sec)  ' + str(rt60_sabine))
        print('Critical distance (m)        ' + str(d_critical))
        print('Mean free path (m)           ' + str(d_mfpath))

    return rt60_sabine, d_critical, d_mfpath
