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
#   @file   compute_echograms.py
#   @author Andrés Pérez-López
#   @date   30/07/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import numpy as np
from .echogram import Echogram
from .image_source_method import ims_coreMtx
from .rec_module import rec_module_mic, rec_module_sh
from .absorption_module import apply_absorption

def compute_echograms_array(room, src, rec, abs_wall, limits):
    """
    Todo
    :param room:
    :param src:
    :param rec:
    :param abs_wall:
    :param limits:
    :return:
    """
    raise NotImplementedError


# TODO: unify abs_wall/alpha arguments (absorption, absorption ratio, diractivity)
def compute_echograms_mic(room, src, rec, abs_wall, limits, mic_specs):
    """
    TODO
    :param room:
    :param src:
    :param rec:
    :param abs_wall:
    :param limits:
    :param mic_specs:
    :return:
    """
    nRec = rec.shape[0]
    nSrc = src.shape[0]
    nBands = abs_wall.shape[0]

    # Limit the RIR by reflection order or by time-limit
    # TODO: expose type as argument
    type = 'maxTime'

    # Compute echogram due to pure propagation (frequency-independent)
    echograms = np.empty((nSrc, nRec), dtype=Echogram)
    for ns in range(nSrc):
        for nr in range(nRec):
            print('Compute echogram: Source ' + str(ns) + ' - Receiver ' + str(nr))
            # Compute echogram
            echograms[ns, nr] = ims_coreMtx(room, src[ns,:], rec[nr,:], type, np.max(limits))

    print('Apply receiver direcitivites')
    rec_echograms = rec_module_mic(echograms, mic_specs)

    abs_echograms = np.empty((nSrc, nRec, nBands), dtype=Echogram)
    # Apply boundary absorption
    for ns in range(nSrc):
        for nr in range(nRec):
            print ('Apply absorption: Source ' + str(ns) + ' - Receiver ' + str(nr))
            # Compute echogram
            abs_echograms[ns, nr, :] = apply_absorption(rec_echograms[ns, nr], abs_wall, limits)

    # TODO: validate return
    # return abs_echograms, rec_echograms, echograms
    return abs_echograms


def compute_echograms_sh(room, src, rec, abs_wall, limits, rec_orders):
    """
    TODO
    :param room:
    :param src:
    :param rec:
    :param abs_wall:
    :param limits:
    :param rec_orders:
    :return:
    """

    nRec = rec.shape[0]
    nSrc = src.shape[0]
    nBands = abs_wall.shape[1]

    # Limit the RIR by reflection order or by time-limit
    type = 'maxTime'
    # Compute echogram due to pure propagation (frequency-independent)
    echograms = np.empty((nSrc, nRec), dtype=Echogram)
    for ns in range(nSrc):
        for nr in range(nRec):
            print('Compute echogram: Source ' + str(ns) + ' - Receiver ' + str(nr))
            # Compute echogram
            echograms[ns, nr] = ims_coreMtx(room, src[ns,:], rec[nr,:], type, np.max(limits))

    print('Apply SH directivites')
    rec_echograms = rec_module_sh(echograms, rec_orders)

    abs_echograms = np.empty((nSrc, nRec, nBands), dtype=Echogram)
    # Apply boundary absorption
    for ns in range(nSrc):
        for nr in range(nRec):
            print ('Apply absorption: Source ' + str(ns) + ' - Receiver ' + str(nr))
            # Compute echogram
            abs_echograms[ns, nr] = apply_absorption(rec_echograms[ns, nr], abs_wall, limits)

    # TODO: validate return
    # return abs_echograms, rec_echograms, echograms
    return abs_echograms

