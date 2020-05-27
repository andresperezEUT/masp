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
#   @file   image_source_method.py
#   @author Andrés Pérez-López
#   @date   30/07/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import numpy as np

from masp.validate_data_types import _validate_ndarray_1D, _validate_int, _validate_number, \
    _validate_echogram, _validate_string
from .echogram import Echogram
from masp.utils import C, c


def ims_coreMtx(room, source, receiver, type, typeValue):
    """
    Compute echogram by image source method.

    Parameters
    ----------
    room : ndarray
        Room dimensions in cartesian coordinates. Dimension = (3) [x, y, z].
    src : ndarray
        Source position in cartesian coordinates. Dimension = (3) [x, y, z].
    rec : ndarray
        Receiver position in cartesian coordinates. Dimension = (3) [x, y, z].
    type : str
        Restriction type: 'maxTime' or 'maxOrder'
    typeValue: int or float
        Value of the chosen restriction.

    Returns
    -------
    reflections : echogram
        An Echogram instance.

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    `src` and `rec` positions are specified from the left ground corner
    of the room, using a left-handed coordinate system.
    `room` refers to the wall dimensions.
    Therefore, their values should be positive and smaller than room dimensions.

              _____    _
             |     |   |
             |     |   |
           x ^     |   | l = r[0]
             |     |   |
             |     |   |
             o---->    -
                  y
             |-----|
                w = r[1]

    """

    _validate_ndarray_1D('room', room, size=C, positive=True)
    _validate_ndarray_1D('source', source, size=C, positive=True, limit=[np.zeros(C),room])
    _validate_ndarray_1D('receiver', receiver, size=C, positive=True, limit=[np.zeros(C),room])
    _validate_string('type', type, choices=['maxTime', 'maxOrder'])

    # Room dimensions
    l, w, h = room

    # Move source origin to the centrer of the room
    src = np.empty(C)
    src[0] = source[0] - l / 2
    src[1] = w / 2 - source[1]
    src[2] = source[2] - h / 2

    # Move receiver origin to the centrer of the room
    rec = np.empty(C)
    rec[0] = receiver[0] - l / 2
    rec[1] = w / 2 - receiver[1]
    rec[2] = receiver[2] - h / 2

    if type is 'maxOrder':
        maxOrder = typeValue
        echogram = ims_coreN(room, src, rec, maxOrder)
    elif type is 'maxTime':
        maxDelay = typeValue
        echogram = ims_coreT(room, src, rec, maxDelay)

    # Sort reflections according to propagation time
    idx = np.argsort(echogram.time)
    echogram.time = echogram.time[idx]
    echogram.value = echogram.value[idx]
    echogram.order = echogram.order[idx, :]
    echogram.coords = echogram.coords[idx, :]

    return echogram


def ims_coreN(room, src, rec, N):
    """
    Compute echogram by image source method, under reflection order restriction

    Parameters
    ----------
    room : ndarray
        Room dimensions in cartesian coordinates. Dimension = (3) [x, y, z].
    src : ndarray
        Source position in cartesian coordinates. Dimension = (3) [x, y, z].
    rec : ndarray
        Receiver position in cartesian coordinates. Dimension = (3) [x, y, z].
    N : int
        Maximum reflection order.

    Returns
    -------
    reflections : echogram
        An Echogram instance.

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    `src` and `rec` positions are specified from a right-handed coordinate system
     placed at the center of the room, with +x facing front, and +y facing left.
     (as opposite to `ims_coreMtx`).
     However, `room` refer to the wall dimensions.
     Therefore, given values must be in the range +-room[i]/2.

                ^x
              __|__    _
             |  |  |   |
             |  |  |   |
          y<----o  |   | l = r[0]
             |     |   |
             |     |   |
             |_____|   -

             |-----|
                w = r[1]

    """

    _validate_ndarray_1D('room', room, size=C, positive=True)
    _validate_ndarray_1D('source', src, size=C, limit=[-room/2,room/2])
    _validate_ndarray_1D('receiver', rec, size=C, limit=[-room/2,room/2])
    _validate_int('N', N, positive=True)

    # i,j,k indices for calculation in x,y,z respectively
    r = np.arange(-N, N+1)
    xx, yy, zz = np.meshgrid(r, r, r)
    # Vectorize (kind of empirical...)
    i = zz.reshape(zz.size)
    j = xx.reshape(xx.size)
    k = yy.reshape(yy.size)
    # Compute total order and select only valid incides up to order N
    s_ord = np.abs(i) + np.abs(j) + np.abs(k)
    i = i[s_ord <= N]
    j = j[s_ord <= N]
    k = k[s_ord <= N]

    # Image source coordinates with respect to receiver
    s_x = i*room[0] + np.power(-1.,i)*src[0] - rec[0]
    s_y = j*room[1] + np.power(-1.,j)*src[1] - rec[1]
    s_z = k*room[2] + np.power(-1.,k)*src[2] - rec[2]

    # Distance
    s_d = np.sqrt(np.power(s_x,2) + np.power(s_y,2) + np.power(s_z,2))
    # Reflection propagation time
    s_t = s_d/c
    # Reflection propagation attenuation - if distance is <1m
    # set at attenuation at 1 to avoid amplification
    s_att = np.empty(s_d.size)
    s_att[s_d <= 1] = 1
    s_att[s_d > 1] = 1./s_d[s_d > 1]

    # Write to echogram structure
    reflections = Echogram(value=s_att[:, np.newaxis],
                           time=s_t,
                           order=np.asarray(np.stack([i, j, k], axis=1), dtype=int),
                           coords=np.stack([s_x, s_y, s_z], axis=1))

    return reflections


def ims_coreT(room, src, rec, maxTime):
    """
    Compute echogram by image source method, under maxTime restriction

    Parameters
    ----------
    room : ndarray
        Room dimensions in cartesian coordinates. Dimension = (3) [x, y, z].
    src : ndarray
        Source position in cartesian coordinates. Dimension = (3) [x, y, z].
    rec : ndarray
        Receiver position in cartesian coordinates. Dimension = (3) [x, y, z].
    maxTime : float
        Maximum echogram computation time.

    Returns
    -------
    reflections : echogram
        An Echogram instance.

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    `src` and `rec` positions are specified from a right-handed coordinate system
     placed at the center of the room, with +x facing front, and +y facing left.
     (as opposite to `ims_coreMtx`).
     However, `room` refer to the wall dimensions.
     Therefore, given values must be in the range +-room[i]/2.

                ^x
              __|__    _
             |  |  |   |
             |  |  |   |
          y<----o  |   | l = r[0]
             |     |   |
             |     |   |
             |_____|   -

             |-----|
                w = r[1]

    """

    _validate_ndarray_1D('room', room, size=C, positive=True)
    _validate_ndarray_1D('source', src, size=C, limit=[-room/2,room/2])
    _validate_ndarray_1D('receiver', rec, size=C, limit=[-room/2,room/2])
    _validate_number('maxTime', maxTime, positive=True)

    # Find order N that corresponds to maximum distance
    d_max = maxTime * c
    Nx = np.ceil(d_max / room[0])
    Ny = np.ceil(d_max / room[1])
    Nz = np.ceil(d_max / room[2])

    # i, j, k indices for calculation in x, y, z respectively
    rx = np.arange(-Nx, Nx + 1)
    ry = np.arange(-Ny, Ny + 1)
    rz = np.arange(-Nz, Nz + 1)
    xx, yy, zz = np.meshgrid(rx, ry, rz)
    # Vectorize (transpose idx due to matlab/python variations on matrix handling)
    i = xx.transpose(2,0,1).flatten()
    j = yy.transpose(2,0,1).flatten()
    k = zz.transpose(2,0,1).flatten()
    # Image source coordinates with respect to receiver
    s_x = i*room[0] + np.power(-1.,i)*src[0] - rec[0]
    s_y = j*room[1] + np.power(-1.,j)*src[1] - rec[1]
    s_z = k*room[2] + np.power(-1.,k)*src[2] - rec[2]
    # Distance
    s_d = np.sqrt(np.power(s_x,2) + np.power(s_y,2) + np.power(s_z,2))

    # Bypass image sources with d > dmax
    i = i[s_d < d_max]
    j = j[s_d < d_max]
    k = k[s_d < d_max]
    s_x = s_x[s_d < d_max]
    s_y = s_y[s_d < d_max]
    s_z = s_z[s_d < d_max]
    s_d = s_d[s_d < d_max]

    # Reflection propagation time
    s_t = s_d/c
    # Reflection propagation attenuation - if distance is <1m
    # set at attenuation at 1 to avoid amplification
    s_att = np.zeros(s_d.size)
    s_att[s_d <= 1] = 1
    s_att[s_d > 1] = 1./s_d[s_d > 1]

    # Write to echogram structure
    reflections = Echogram(value=s_att[:, np.newaxis],
                           time=s_t,
                           order=np.asarray(np.stack([i, j, k], axis=1), dtype=int),
                           coords=np.stack([s_x, s_y, s_z], axis=1))

    return reflections
