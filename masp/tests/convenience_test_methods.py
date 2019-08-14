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
#   run it with `python3 -m pytest --cov-report term-missing --cov=masp -v -s masp/tests/`
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
# matlab session will automatically close when shoebox_room_sim properly finish,
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

# Path to tmp folder for Matlab struct storage and evaluation
tmp_path = os.path.abspath("./masp/tests/tmp")


def get_parameters(params, t):
    """
    TODO
    get a params dictionary, with each key a parameter name,
    and value a list of values of lenght T.
    :param params:
    :param t:
    :return:
    """
    num_params = len(params)
    print('')
    print('-----------------------------------------------')
    print('                  t='+str(t))
    p = []
    for p_idx in range(num_params):
        dict_tuple = list(params.items())[p_idx]
        key = dict_tuple[0]
        value = dict_tuple[1]
        p.append(value[t])
        print(key, value[t])
    print('-----------------------------------------------')
    return p


def generate_random_echogram():
    room = np.random.random(C) * 5 + 5
    src = np.random.random(C) * 5 - 2.5
    rec = np.random.random(C) * 5 - 2.5
    N = np.random.randint(20)
    echo = masp.srs.ims_coreN(room, src, rec, N)
    _validate_echogram(echo)
    return echo

def generate_random_echogram_array(nSrc, nRec):
    echogram_array = np.empty((nSrc, nRec), dtype=masp.srs.Echogram)
    for src in range(nSrc):
        for rec in range(nRec):
            echogram_array[src, rec] = generate_random_echogram()
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

def compare_echograms(np_res, ml_res):
    if not np.allclose(np.asarray(ml_res['time']).squeeze(), np_res.time.squeeze()): raise_error()
    if not np.allclose(np.asarray(ml_res['value']).squeeze(), np_res.value): raise_error()
    if not np.allclose(ml_res['order'], np_res.order): raise_error()
    if not np.allclose(ml_res['coords'], np_res.coords): raise_error()

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
        if not np.allclose(ml_res['value'][idx].squeeze(), np_res[idx].value.squeeze()): raise_error()
        if not np.allclose(ml_res['order'][idx], np_res[idx].order): raise_error()
        if not np.allclose(ml_res['coords'][idx], np_res[idx].coords): raise_error()


def numeric_assert(ml_method, np_method, *args, nargout=0, write_file=False, namespace=None):

    # Convert arguments to required data types
    ml_args = []
    np_args = []
    for arg in args:

        if isinstance(arg, list):
            ml_args.append(matlab.double(arg))
            np_args.append(np.array(arg))

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

        # Array of echograms: write to file and append path to read
        elif isinstance(arg, np.ndarray):
            # Parse array and convert to array of dicts
            s0, s1 = arg.shape
            ml_echogram_array = np.empty((s0, s1), dtype=masp.srs.Echogram)
            for i in range(s0):
                for j in range(s1):
                    echo = arg[i,j]
                    ml_echo = {}
                    ml_echo['time'] = echo.time[:, np.newaxis].tolist()
                    ml_echo['value'] = echo.value[:, np.newaxis].tolist()
                    ml_echo['order'] = echo.order.tolist()
                    ml_echo['coords'] = echo.coords.tolist()
                    ml_echogram_array[i, j] = ml_echo

            # Save resulting array to tmp folder, and add the path to the method args
            ml_echogram_array_path = os.path.join(tmp_path, ml_method+'_arg.mat')
            if os.path.isfile(ml_echogram_array_path):
                os.remove(ml_echogram_array_path)
            savemat(ml_echogram_array_path, {'arg': ml_echogram_array})
            ml_args.append(ml_echogram_array_path)

            # For python the argument is fine
            np_args.append(arg)

        elif isinstance(arg, str):
            ml_args.append(arg)
            np_args.append(arg)

        elif arg is None:
            pass

        else:
            print(type(arg))
            raise NotImplementedError

    # Run method
    if namespace is None:
        np_res = getattr(masp, np_method)(*np_args)
    else:
        np_res = getattr(getattr(masp, namespace), np_method)(*np_args)

    if write_file:
        ml_args.insert(0, tmp_path)
        ml_method_test = ml_method+'_test'
        getattr(eng, ml_method_test)(*ml_args, nargout=nargout)
    else:
        ml_res = getattr(eng, ml_method)(*ml_args, nargout=nargout)


    # Write check: matlab-python interface is not able to pass struct arrays.
    # Therefore, we call the `test` version of the matlab function,
    # which writes the resulting data as .mat v6 into our tmp folder,
    # with the same name as the matlab function called
    # Then, we can open it using scipy.io
    if write_file:
        # TODO: WHEN MULTIPLE OUTPUT FILES, NAME THEM SIMILARLY AND FIND THEM BY NAME HERE
        tmp_file_path = os.path.join(tmp_path, ml_method+'_test.mat')
        mat_file = loadmat(tmp_file_path)

        # mat_file is a dictionary containing the data along some metadata...
        # The actual data is stored in the only "non-private" key...
        ml_res = next( mat_file[k] for k in mat_file.keys() if not k.startswith("__") )

        # Actually compare them
        compare_echogram_arrays(np_res, ml_res)

        # Remove used file
        os.remove(tmp_file_path)
        # Also remove matlab arg file, in case
        if 'ml_echogram_array_path' in locals():
            if os.path.isfile(ml_echogram_array_path):
                os.remove(ml_echogram_array_path)

    # Regular check: matlab result data is saved into ml_res
    else:
        if nargout == 1:
            if isinstance(np_res, np.ndarray):
                if not np.allclose(ml_res, np_res): raise_error()

            elif isinstance(np_res, float):
                if ml_res != np_res: raise_error()

            elif isinstance(np_res, masp.srs.Echogram):
                compare_echograms(np_res, ml_res)

            else:
                print(type(np_res))
                raise_error(NotImplementedError)

        else:
            for arg_idx in range(nargout):

                if isinstance(np_res[arg_idx], np.ndarray):
                    if not np.allclose(ml_res[arg_idx], np_res[arg_idx]): raise_error()

                elif isinstance(np_res[arg_idx], float):
                    if ml_res[arg_idx] != np_res[arg_idx]: raise_error()

                elif isinstance(arg, masp.srs.Echogram):
                    compare_echograms(np_res, ml_res)

                else:
                    print(type(np_res))
                    raise_error(NotImplementedError)

