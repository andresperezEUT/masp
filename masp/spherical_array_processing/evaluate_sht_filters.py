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
#   @file   evaluate_sht_filters.py
#   @author Andrés Pérez-López
#   @date   01/10/2019
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import numpy as np
import matplotlib.pyplot as plt
from masp.validate_data_types import _validate_int, _validate_ndarray_1D, \
    _validate_ndarray_2D, _validate_ndarray_3D, _validate_boolean, _validate_number


def evaluate_sht_filters(M_mic2sh, H_array, fs, Y_grid, w_grid=None, plot=False):
    """
    Evaluate frequency-dependent performance of SHT filters.

    Parameters
    ----------
    M_mic2sh : ndarray
        SHT filtering matrix produced by one of the methods included in the library.
        Dimension = ( (order+1)^2, nMics, nBins ).
    H_array : ndarray, dtype = 'complex'
         Modeled or measured spherical array responses in a dense grid of `nGrid` directions.
        Dimension = ( nBins, nMics, nGrid ).
    fs : int
        Target sampling rate.
    Y_grid : ndarray
        Spherical harmonics matrix for the `nGrid` directions of the evaluation grid.
        Dimension = ( nGrid, (order+1)^2 ).
    w_grid : ndarray, optional
        Vector of integration weights for the grid points.
        Dimension = ( nGrid ).
    plot : bool, optional
        Plot responses. Default to false.

    Returns
    -------
    cSH : ndarray, dtype = 'complex'
        Spatial correlation coefficient, for each SHT order and frequency bin.
        Dimension = ( nBins, order+1 ).
    lSH : ndarray
        Level difference, for each SHT order, for each SHT order and frequency bin.
        Dimension = ( nBins, order+1 ).
    WNG : ndarray
        Maximum amplification of all output SH components.
        Dimension = ( nBins ).

    Raises
    -----
    TypeError, ValueError: if method arguments mismatch in type, dimension or value.

    Notes
    -----
    The SHT filters can be evaluated in terms of how ideal are the SH
    components that they generate. The evaluation here follows the metrics
    introduced in

        Moreau, S., Daniel, J., Bertet, S., 2006,
        `3D sound field recording with higher order ambisonics-objectiv
        measurements and validation of spherical microphone.`
        In Audio Engineering Society Convention 120.

    These are a) the spatial correlation coefficient between each ideal
    spherical harmonic and the reconstructed pattern, evaluated at a dense
    grid of directions, b) level difference between the mean spatial power
    of the reconstructed pattern (diffuse power) over the one from an ideal
    SH component. Ideally, correlaiton should be close to one, and the
    level difference should be close to 0dB.

    Additionally, the maximum amplification of all output SH components is
    evaluated, through the maximum eigenvalue of the filtering matrix.

    Due to the matrix nature of computations,
    the minimum valid value for `nMics` and `nGrid` is 2.
    """

    _validate_ndarray_3D('M_mic2sh', M_mic2sh)
    n_sh = M_mic2sh.shape[0]
    order_sht = int(np.sqrt(n_sh) - 1)
    nMics = M_mic2sh.shape[1]
    _validate_number('nMics', nMics, limit=[2, np.inf])
    nBins = M_mic2sh.shape[2]

    _validate_ndarray_3D('H_array', H_array, shape0=nBins, shape1=nMics)
    nGrid = H_array.shape[2]
    _validate_number('nGrid', nGrid, limit=[2, np.inf])

    _validate_ndarray_2D('Y_grid', Y_grid, shape0=nGrid, shape1=n_sh)

    if w_grid is None:
        w_grid = 1/nGrid*np.ones(nGrid)
    _validate_ndarray_1D('w_grid', w_grid, size=nGrid)

    _validate_int('fs', fs, positive=True)
    if plot is not None:
        _validate_boolean('plot', plot)

    nFFT = 2 * (nBins - 1)
    f = np.arange(nFFT // 2 + 1) * fs / nFFT

    # Compute spatial correlations and integrated level difference between
    # ideal and reconstructed harmonics
    cSH = np.empty((nBins, order_sht+1), dtype='complex')
    lSH = np.empty((nBins, order_sht+1))
    # rSH = np.empty((nBins, order_sht+1))
    for kk in range(nBins):
        H_kk = H_array[kk,:,:]
        y_recon_kk = np.matmul(M_mic2sh[:,:, kk], H_kk)
        for n in range(order_sht+1):
            cSH_n = 0  # spatial correlation (mean per order)
            lSH_n = 0  # diffuse level difference (mean per order)
            # rSH_n = 0  # mean level difference (mean per order)
            for m in range(-n, n+1):
                q = np.power(n, 2) + n + m
                y_recon_nm = y_recon_kk[q,:].T
                y_ideal_nm = Y_grid[:, q]
                cSH_nm = np.matmul((y_recon_nm * w_grid).conj(), y_ideal_nm) / np.sqrt( np.matmul((y_recon_nm*w_grid).conj(), y_recon_nm ))
                cSH_n = cSH_n + cSH_nm
                lSH_nm = np.real(np.matmul((y_recon_nm * w_grid).conj(), y_recon_nm ))
                lSH_n = lSH_n + lSH_nm
                # rSH_nm = np.sum(np.power(np.abs(y_recon_nm - y_ideal_nm), 2) * w_grid)
                # rSH_n = rSH_n + rSH_nm;
            cSH[kk, n] = cSH_n / (2 * n + 1)
            lSH[kk, n] = lSH_n / (2 * n + 1)
            # rSH[kk, n] = rSH_n / (2 * n + 1)

    # Maximum noise amplification of all filters in matrix
    WNG = np.empty(nBins)
    for kk in range(nBins):
        # TODO: Matlab implementation warns when M matrix is complex, e.g. TEST_SCRIPTS l. 191-199
        # Avoid ComplexWarning: imaginary parts appearing due to numerical precission
        eigM = np.real(np.linalg.eigvals(np.matmul(M_mic2sh[:,:,kk].T.conj(), M_mic2sh[:,:,kk])))
        WNG[kk] = np.max(eigM)

    # Plots
    if plot:
        str_legend = [None]*(order_sht+1)
        for n in range(order_sht+1):
            str_legend[n] = str(n)

        plt.figure()
        plt.subplot(311)
        plt.semilogx(f, np.abs(cSH))
        plt.grid()
        plt.legend(str_legend)
        plt.axis([50, 20000, 0, 1])
        plt.title('Spatial correlation')

        plt.subplot(312)
        plt.semilogx(f, 10 * np.log10(lSH))
        plt.grid()
        plt.legend(str_legend)
        plt.axis([50, 20000, -30, 10])
        plt.title('Level correlation')

        plt.subplot(313)
        plt.semilogx(f, 10 * np.log10(WNG))
        plt.grid()
        plt.xlim([50, 20000])
        plt.title('Maximum amplification')
        plt.xlabel('Frequency (Hz)')

        # plt.subplot(414)
        # plt.semilogx(f, 10 * np.log10(rSH))
        # plt.grid()
        # plt.xlim([50, 20000])
        # plt.title('MSE')
        # plt.xlabel('Frequency (Hz)')

        plt.show()

    return cSH, lSH, WNG
