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
#   @file   plot_functions.py
#   @author Andrés Pérez-López
#   @date   27/09/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from masp.utils import sph2cart, C
from masp.validate_data_types import _validate_ndarray_2D


def plot_mic_array(capsule_positions):
    """
    Depict a 3D plot of a given microphone array geometry.

    Parameters
    ----------
    capsule_positions : ndarray
      Spherical coordinates (in radians) of the capsules. Dimension = (nMic, C)

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.



    """

    _validate_ndarray_2D('capsule_positions', capsule_positions, shape1=C)

    nMic = capsule_positions.shape[0]
    R = capsule_positions[-0,-1]
    cart = sph2cart(capsule_positions)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # sphere code adapted from https://stackoverflow.com/questions/33551103/plotting-a-sphere-mesh-with-matplotlib
    u = np.linspace(0, 2 * np.pi, 30)
    v = np.linspace(0, np.pi, 15)
    x = 0.99*R * np.outer(np.cos(u), np.sin(v))
    y = 0.99*R * np.outer(np.sin(u), np.sin(v))
    z = 0.99*R * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z, rstride=1, cstride=1, color='w', shade=1, alpha=0.5)

    # TODO: draw axes
    # ax.scatter(cart[:, 0], cart[:, 1], cart[:, 2], s=100)
    # Ffor each given sphere, shift the scaled unit sphere by the
    # location of the sphere and plot

    # Draw points and capsule numbers
    for nm in range(nMic):
        ax.scatter(cart[nm, 0], cart[nm, 1], cart[nm, 2], marker='o')
        ax.text(cart[nm, 0], cart[nm, 1], cart[nm, 2], s=str(nm+1))
