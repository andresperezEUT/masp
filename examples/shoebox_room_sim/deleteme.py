from masp import shoebox_room_sim as srs
import numpy as np
import time

# room = np.array([ 4, 7, 10])
# src = np.array([2.1571,   -2.2931,   -2.2451])
# rec = np.array([1.5632,   -0.3969,    1.6549])
# maxTime =  0.4047
#
# a = time.time()
# r = srs.ims_coreT(room, src, rec, maxTime)
# b = time.time()
# print(b-a)
# print(np.size(r.value))


import masp
import os.path
from scipy.io import loadmat, savemat
from numpy.core.records import fromarrays


tmp_path = os.path.abspath("/Users/andres.perez/source/masp/masp/tests/tmp")


### 1x1 STRUCT

# generate_random_echogram():
room = np.random.random(3) * 5 + 5
src = np.random.random(3) * 5 - 2.5
rec = np.random.random(3) * 5 - 2.5
N = np.random.randint(20)
echo = masp.srs.ims_coreN(room, src, rec, N)


## METHOD 1
## Note: this approach is working as well!
# # Convert 1d python arrays to column vectors in python
# echo.time = echo.time[:, np.newaxis]
# echo.value = echo.value[:, np.newaxis]
#
# # This is properly loaded in matlab as a 1x1 struct!
# savemat(os.path.join(tmp_path, 'echo1x1.mat'), {'echo': echo})

## METHOD 2
# Somehow the above approach is working, but the proper approach
# would be to have echo directly as a dict
echo_dict = {}
echo_dict['time'] = echo.time[:, np.newaxis]
echo_dict['value'] = echo.value[:, np.newaxis]
echo_dict['order'] = echo.order
echo_dict['coords'] = echo.coords

# This is properly loaded in matlab as a 1x1 struct!
# Note: saving echo_dict directly will cause matlab
# to expand the individual variables (instead of creating a struct)
savemat(os.path.join(tmp_path, 'echo1x1.mat'), {'echo': echo_dict})



## Load the same 1x1 struct
# loadmat yields a dictionary with some stuff, but the interesting
# data is given in the value corresponding to the variable name's key
x = loadmat(os.path.join(tmp_path, 'echo1x1.mat'))['echo']

# now x is a structured nd rray (https://docs.scipy.org/doc/numpy/user/basics.rec.html)
# x.dtype -> dtype([('value', 'O'), ('time', 'O'), ('order', 'O'), ('coords', 'O')])

# furthermore, each of the 'fields' of x has a shape of (1,1)
assert x['value'].shape == x['coords'].shape == (1,1)
# which suggests a way to handle matlab struct arrays.

# The actual data is stored here at the first element (in 2D)...
print(x['value'][0,0].shape) # Should be (n, 1)
print(x['coords'][0,0].shape) # Should be (n, 3)

assert np.allclose(x['value'][0,0], echo_dict['value'])
assert np.allclose(x['time'][0,0], echo_dict['time'])
assert np.allclose(x['order'][0,0], echo_dict['order'])
assert np.allclose(x['coords'][0,0], echo_dict['coords'])




### NxN STRUCT
# In this case, we create a NxM ndarray of type Object (Echogram in this case),
# and write it directly with `savemat`.
# Matlab will interpret it as a NxM cell, not struct, but the transformation is straightforward.
# This procedure could be standardized for the general case.

def generate_random_echogram():
    room = np.random.random(3) * 5 + 5
    src = np.random.random(3) * 5 - 2.5
    rec = np.random.random(3) * 5 - 2.5
    # N = np.random.randint(20)
    N = 1
    echo = masp.srs.ims_coreN(room, src, rec, N)

    ml_echo = {}
    ml_echo['time'] = echo.time[:, np.newaxis]
    ml_echo['value'] = echo.value[:, np.newaxis]
    ml_echo['order'] = echo.order
    ml_echo['coords'] = echo.coords

    # return echo
    return ml_echo

def generate_random_echogram_array(nSrc, nRec):
    echogram_array = np.empty((nSrc, nRec), dtype=masp.srs.Echogram)
    for src in range(nSrc):
        for rec in range(nRec):
            echogram_array[src, rec] = generate_random_echogram()
    return echogram_array


nSrc = np.random.randint(1, 10)
nRec = np.random.randint(1, 10)
struct_array = generate_random_echogram_array(nSrc, nRec)
savemat(os.path.join(tmp_path,'struct_array.mat'), {'struct_array': struct_array})





