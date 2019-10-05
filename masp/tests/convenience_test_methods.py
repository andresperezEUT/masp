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
#   @file   convenience_test_methods.py
#   @author Andrés Pérez-López
#   @date   31/07/2019
#
#   run it with `python3 -m pytest --cov-report term-missing --cov=masp -v -s -x masp/tests/`
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import matlab.engine
import masp
import numpy as np
import os.path
from scipy.io import loadmat, savemat
from masp.validate_data_types import _validate_echogram, _validate_echogram_array
from masp.utils import C


# Start matlab
# matlab session will automatically close when test session properly finishes,
# or alternatively on AssertionError handling

eng = matlab.engine.start_matlab()

# Path to matlab code to test
matlab_path = '/Users/andres.perez/source/MATLAB/polarch/'
eng.addpath(matlab_path)
# Shoebox library
eng.addpath(os.path.join(matlab_path,'shoebox-roomsim-master'))
eng.addpath(os.path.join(matlab_path,'shoebox-roomsim-master', 'test'))

# Spherical Harmonic Transform library
eng.addpath(os.path.join(matlab_path,'Spherical-Harmonic-Transform-master'))
# Array Response Simulator library
eng.addpath(os.path.join(matlab_path,'Array-Response-Simulator-master'))
# Spherical Array Processing library
eng.addpath(os.path.join(matlab_path,'Spherical-Array-Processing-master'))

# Path to tmp folder for Matlab struct storage and evaluation
tmp_path = os.path.abspath("./masp/tests/tmp")

# Define echogram dtypes...
echogram_dtype = np.dtype([('time', 'O'), ('value', 'O'), ('order', 'O'), ('coords', 'O')])
quantised_echogram_dtype = np.dtype([('value', 'O'), ('time', 'O'), ('isActive', 'O')])

def get_parameters(params, t, verbose=True):
    """
    Convenience method for the numeric testing/validation system.

    Params is a dictionary, each key being a parameter name,
    and each value a list of values of length T (the number of trials).

    The method returns a list of the argument values for test number `t` < T.
    """

    def print_v(str):
        if verbose:
            print(str)

    num_params = len(params)
    print_v('')
    print_v('-----------------------------------------------')
    print_v('                  t='+str(t))
    p = []
    for p_idx in range(num_params):
        dict_tuple = list(params.items())[p_idx]
        key = dict_tuple[0]
        value = dict_tuple[1]
        p.append(value[t])
        if verbose:
            print(key, value[t])
    print_v('-----------------------------------------------')
    return p


def generate_random_echogram():
    room = np.random.random(C) * 5 + 5
    src = np.random.random(C) * 2.5
    rec = np.random.random(C) * 2.5
    N = np.random.randint(20)+3
    echo = masp.srs.ims_coreMtx(room, src, rec, 'maxOrder', N)
    return echo

def generate_random_echogram_sh(nSH=4):
    room = np.random.random(C) * 5 + 5
    src = np.random.random(C) * 2.5
    rec = np.random.random(C) * 2.5
    N = np.random.randint(20)+3
    echo = masp.srs.ims_coreMtx(room, src, rec, 'maxOrder', N)
    # Expand the value dimension simulating sh echograms
    echo.value = echo.value*np.random.rand(nSH)
    return echo

def generate_random_echogram_array(nSrc, nRec, nBands=None):
    if nBands:
        echogram_array = np.empty((nSrc, nRec, nBands), dtype=masp.srs.Echogram)
    else:
        echogram_array = np.empty((nSrc, nRec), dtype=masp.srs.Echogram)

    for idx in np.ndindex(echogram_array.shape):
        echogram_array[idx] = generate_random_echogram()

    _validate_echogram_array(echogram_array)
    return echogram_array

def generate_random_echogram_array_sh(nSrc, nRec, nBands=None):
    if nBands:
        echogram_array = np.empty((nSrc, nRec, nBands), dtype=masp.srs.Echogram)
    else:
        echogram_array = np.empty((nSrc, nRec), dtype=masp.srs.Echogram)

    nSH = np.random.randint(1,10)
    for idx in np.ndindex(echogram_array.shape):
        echogram_array[idx] = generate_random_echogram_sh(nSH)

    _validate_echogram_array(echogram_array)
    return echogram_array


def generate_random_unit_vector():
    v = np.random.rand(C)
    return (v/np.linalg.norm(v)).tolist()

def generate_random_mic_specs(nRec):
    mic_specs = []
    for rec in range(nRec):
        mic_specs.append(generate_random_unit_vector())
        mic_specs[-1].append(np.random.rand())
    return mic_specs



def raise_error(e=None):
    eng.quit()
    if e is None:
        raise AssertionError
    else:
        raise e

def validate_result(ml_res, np_res, rtol, atol):
    if isinstance(np_res, np.ndarray):
        m = np.asarray(ml_res).squeeze()
        n = np_res.squeeze()
        if not m.shape == n.shape: raise_error()
        if not np.allclose(m, n, rtol=rtol, atol=atol): raise_error()

    elif isinstance(np_res, list):
        # Matlab cell -> python list
        # List elements should be ndarrays
        if not len(np_res) == len(ml_res): raise_error()
        for i in range(len(np_res)):
            m = np.asarray(ml_res[i]).squeeze()
            n = np_res[i].squeeze()
            if not m.shape == n.shape: raise_error()
            if not np.allclose(m, n, rtol=rtol, atol=atol): raise_error()

    elif isinstance(np_res, float):
        if ml_res != np_res: raise_error()

    elif isinstance(np_res, masp.srs.Echogram):
        compare_echograms(np_res, ml_res)

    elif isinstance(np_res, masp.srs.QuantisedEchogram):
        compare_quantised_echograms(np_res, ml_res)

    else:
        print(type(np_res))
        raise_error(NotImplementedError)

def compare_echograms(np_res, ml_res):
    # In python, 'time' is 1D, while others are 2D. Therefore we must use squeeze for comparison
    if not np.allclose(np.asarray(ml_res['time']).squeeze(), np_res.time.squeeze()): raise_error()
    if not np.allclose(ml_res['value'], np_res.value): raise_error()
    if not np.allclose(ml_res['order'], np_res.order): raise_error()
    if not np.allclose(ml_res['coords'], np_res.coords): raise_error()

def compare_quantised_echograms(np_res, ml_res):
    # In python, 'time' is 1D, while others are 2D. Therefore we must use squeeze for comparison
    if not np.allclose(np.asarray(ml_res['time']).squeeze(), np_res.time.squeeze()): raise_error()
    if not np.allclose(ml_res['value'], np_res.value): raise_error()
    if not np.allclose(ml_res['isActive'], np_res.isActive): raise_error()

def compare_echogram_arrays(np_res, ml_res):
    # Compare shapes
    # In the case of `absorption_module`, the resulting matlab abs_echogram
    # might be 2D instead of 3D if nBands==1 (due to struct expansion/filling).
    # So, in this case, just adapt the numpy array for the comparison
    if np_res.ndim == 3 and np_res.shape[-1] == 1:
        np_res = np_res.squeeze(axis=-1)
    if not ml_res.shape == np_res.shape: raise_error()
    # Compare values
    for idx in np.ndindex(np_res.shape):    # multidimensional iterator, valid for 2D and 3D
        if not np.allclose(ml_res['time'][idx].squeeze(), np_res[idx].time.squeeze()): raise_error()
        # if not np.allclose(ml_res['value'][idx].squeeze(), np_res[idx].value.squeeze()): raise_error()
        if not np.allclose(ml_res['value'][idx], np_res[idx].value): raise_error()
        if not np.allclose(ml_res['order'][idx], np_res[idx].order): raise_error()
        if not np.allclose(ml_res['coords'][idx], np_res[idx].coords): raise_error()

def compare_quantised_echogram_arrays(np_res, ml_res):
    # For the moment, this function compares the output of `quantise_echogram()`
    # Compare shapes
    ml_res = ml_res.squeeze()
    np_res = np_res.squeeze()
    if not ml_res.shape == np_res.shape: raise_error()
    # Compare values
    for idx in np.ndindex(np_res.shape):    # multidimensional iterator, valid for 2D and 3D
        if not np.allclose(ml_res['time'][idx].squeeze(), np_res[idx].time.squeeze()): raise_error()
        # if not np.allclose(ml_res['value'][idx].squeeze(), np_res[idx].value.squeeze()): raise_error()
        if not np.allclose(ml_res['value'][idx], np_res[idx].value): raise_error()
        if ml_res['isActive'][idx] != np_res[idx].isActive: raise_error()

def numeric_assert(ml_method, np_method, *args, nargout=0, write_file=False, namespace=None, rtol=1e-5, atol=1e-8):

    # PREPROCESS ARGS --------------------------------------------------------

    # Convert arguments to required data types
    ml_args = []
    np_args = []
    ml_path_arg = False  # flag, active if passing class objects to matlab
    for arg in args:

        if isinstance(arg, list):
            # cast lists to matlab doubles (matrices)
            if np.asarray(arg).dtype != np.dtype('O'):

                # Check if complex, comparing with first element (why not?)
                # TODO: implement in other data types, if needed...
                complex = True if type(np.asarray(arg).flatten()[0]) == np.complex128 else False

                # 1D arrays to matlab column
                if np.asarray(arg).ndim == 1:
                    ml_args.append( matlab.double((np.asarray(arg)[:,np.newaxis]).tolist(), is_complex=complex) )
                else:
                    ml_args.append(matlab.double(arg, is_complex=complex))
                np_args.append(np.array(arg))

        elif isinstance(arg, tuple):
            # automatic tuple->cell conversion in matlab API
            # but still need to convert inner lists to matlab doubles
            list_arg = list(arg)
            for l in range(len(list_arg)):
                list_arg[l] = matlab.double(list_arg[l])
            ml_args.append(tuple(list_arg))

            # in python, that should be a list in the top dimension, and ndarray in the rest
            list_arg = list(arg)
            for l in range(len(list_arg)):
                list_arg[l] = np.asarray(list_arg[l])
            np_args.append(list_arg)

        elif isinstance(arg, (int, float)):
            ml_args.append(np.float(arg))
            np_args.append(arg)

        elif isinstance(arg, masp.srs.Echogram):
            ml_echo = {}
            ml_echo['time'] = matlab.double(arg.time[:,np.newaxis].tolist())
            ml_echo['value'] = matlab.double(arg.value[:,np.newaxis].tolist())
            ml_echo['order'] = matlab.double(arg.order.tolist())
            ml_echo['coords'] = matlab.double(arg.coords.tolist())
            ml_args.append(ml_echo)
            np_args.append(arg)

        elif isinstance(arg, masp.srs.QuantisedEchogram):
            ml_echo = {}
            ml_echo['time'] = matlab.double(arg.time[:,np.newaxis].tolist())
            ml_echo['value'] = matlab.double(arg.value[:,np.newaxis].tolist())
            ml_echo['isActive'] = matlab.double(arg.isActive.tolist())
            ml_args.append(ml_echo)
            np_args.append(arg)

        # Array of echograms: write to file and append path to read
        elif isinstance(arg, np.ndarray) and arg.dtype == np.dtype('O'):
            # Parse array and convert to array of dicts
            if type(arg.flat[0]) == masp.srs.Echogram:  # Check array type
                ml_echogram_array = np.empty(arg.shape, dtype=masp.srs.Echogram)
                for idx in np.ndindex(arg.shape):
                    echo = arg[idx]
                    ml_echo = {}
                    ml_echo['time'] = echo.time[:, np.newaxis].tolist()  # Time is the only 1D field in python
                    ml_echo['value'] = echo.value.tolist()
                    ml_echo['order'] = echo.order.tolist()
                    ml_echo['coords'] = echo.coords.tolist()
                    ml_echogram_array[idx] = ml_echo
            elif type(arg.flat[0]) == masp.srs.QuantisedEchogram:  # Check array type
                ml_echogram_array = np.empty(arg.shape, dtype=masp.srs.QuantisedEchogram)
                for idx in np.ndindex(arg.shape):
                    echo = arg[idx]
                    ml_echo = {}
                    ml_echo['time'] = echo.time[:, np.newaxis].tolist()  # Time is the only 1D field in python
                    ml_echo['value'] = echo.value.tolist()
                    ml_echo['isActive'] = echo.isActive
                    ml_echogram_array[idx] = ml_echo
            else:
                raise TypeError

            # Save resulting array to tmp folder, and add the path to the method args
            ml_echogram_array_path = os.path.join(tmp_path, ml_method+'_arg.mat')
            if os.path.isfile(ml_echogram_array_path):
                os.remove(ml_echogram_array_path)
            savemat(ml_echogram_array_path, {'arg': ml_echogram_array})
            ml_args.append(ml_echogram_array_path)

            # For python the argument is fine
            np_args.append(arg)
            ml_path_arg = True

        elif isinstance(arg, str):
            ml_args.append(arg)
            np_args.append(arg)

        elif arg is None:
            pass

        else:
            print(type(arg))
            raise NotImplementedError


    # RUN METHODS --------------------------------------------------------

    # # Run python method
    # if namespace is None:
    #     np_res = getattr(masp, np_method)(*np_args)
    # else:
    #     np_res = getattr(getattr(masp, namespace), np_method)(*np_args)

    # Run matlab method.
    # Specific matlab test file exists when :
    # 1). passing a file path to be loaded inside matlab (`ml_path_arg=True`)
    # 2). returning a file path from matlab to be read in python (`write_file=True`)

    if write_file:
        ml_args.insert(0, tmp_path)  # Add first argument: path to save output file
        ml_method_test = ml_method + '_test'
        getattr(eng, ml_method_test)(*ml_args, nargout=nargout)
    else:
        if ml_path_arg:
            ml_method += '_test'
        ml_res = getattr(eng, ml_method)(*ml_args, nargout=nargout)

    # Run python method
    if namespace is None:
        np_res = getattr(masp, np_method)(*np_args)
    else:
        np_res = getattr(getattr(masp, namespace), np_method)(*np_args)


    # VALIDATE OUTPUT --------------------------------------------------------

    # Write check: matlab-python interface is not able to pass struct arrays.
    # Therefore, we call the `test` version of the matlab function,
    # which writes the resulting data as .mat v6 into our tmp folder,
    # with the same name as the matlab function called
    # Then, we can open it using scipy.io
    if write_file:
        tmp_file_path = os.path.join(tmp_path, ml_method+'_test.mat')
        mat_file = loadmat(tmp_file_path)

        # mat_file is a dictionary containing the data along some metadata...
        # The actual data is stored in the only "non-private" key...
        ml_res = next( mat_file[k] for k in mat_file.keys() if not k.startswith("__") )

        # Actually compare them
        # First check what kind of result we have...
        if ml_res.dtype == echogram_dtype:
            compare_echogram_arrays(np_res, ml_res)
        elif ml_res.dtype == quantised_echogram_dtype:
            compare_quantised_echogram_arrays(np_res, ml_res)
        else:
            raise TypeError

        # Remove used tmp file
        os.remove(tmp_file_path)

    # Regular check: matlab result data is saved into ml_res
    else:
        if nargout == 1:
            validate_result(ml_res, np_res, rtol, atol)
        else:
            for arg_idx in range(nargout):
                validate_result(ml_res[arg_idx], np_res[arg_idx], rtol, atol)

    # Remove matlab arg file, in case
    if 'ml_echogram_array_path' in locals():
        if os.path.isfile(ml_echogram_array_path):
            os.remove(ml_echogram_array_path)


