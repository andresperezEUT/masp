#MASP ::: Multichannel Acoustic Signal Processing Library for Python

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![PyPI](https://img.shields.io/badge/python-3.6-blue.svg)]()


Python port of the following wonderful [MATLAB libraries from A. Politis](https://github.com/polarch/)>

* _shoebox-roomsim_: Fast Implementation of the Image Source Method 
    * Convex 3D rooms.
    * Arbitrary number of sources and receivers, with arbitrary positions, orientations and directivities.
    * ISM expansion limited by order or time.
    * Frequency-dependent wall absorption.
    * RIR spherical harmonic expansion.

* _array-response-simulator_: Simulation of spherical microphone arrays
    * Rigid/open configurations.
    * Scattering simulation.
    * Arbitrary capsule distances, positions and directivities.

* _spherical-array-processing_: Methods for the transformation and analysis of signals measured with a spherical microphone array
    * A2B conversion with theoretical or measured filters.
    * Signal-independent beamforming. (TODO)
    * Signal-dependent and adaptive beamforming. (TODO)
    * Direction of Arrival estimation. (TODO)
    * Diffuseness estimation. (TODO)
    
* _spherical-harmonic-transform_: Mathematical convenience tools.

Tested in OSX python3.7. 
It should work on python versions >=3.6. 


## Installation

Using pypi:
`pip install masp`


## Examples

Example implementations can be found in the `/examples` folder. 

## Test
Test files are located in `masp/tests/`.
Numeric assertion of the algorithms is performed using the python-Matlab wrapper.
Each method is compared against the Matlab reference implementation. 

Therefore, in order to run the tests, a working copy of Matlab is required (with the signal processing toolbox), along with the [Matlab engine API for Python](https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html).

__________________

### Comparison with pyroomacoustics


| Feature 	                                        | masp | pra |	
|---	                                            |:---: |:---:|
|                                                   |      |  	 | 
| masp.shoebox_room_model                           |      |  	 | 
| --------------------------------------------------|      |  	 | 
| Convex 3D room 	                                | ☑    | ☑ 	 | 
| Non-convex 3D room                                | -    | ☑ 	 | 
| Arbitrary #sources                                | ☑    | ☑ 	 | 
| Arbitrary #receivers, arrays                      | ☑    | ☑ 	 | 
| ISM by max_order                                  | ☑    | ☑ 	 | 
| ISM by max_time                                   | ☑    | - 	 | 
| Wall absorption                                   | ☑    | ☑ 	 | 
| Frequency-dependent absorption                    | ☑    | - 	 | 
| Plot methods                                      | -    | ☑ 	 | 
| RIR rendering                                     | ☑    | ☑ 	 | 
| Audio simulation                                  | ☑    | ☑ 	 | 
| Acoustic descriptor estimation                    | ☑    | - 	 | 
| Microphone orientation                            | ☑    | - 	 | 
| Custom microphone directivity                     | ☑    | - 	 | 
| RIR Spherical Harmonic Expansion                  | ☑    | - 	 | 
|                                                   |      |  	 | 
|                                                   |      |  	 | 
| masp.array_response_simulator                     |      |  	 | 
| --------------------------------------------------|      |  	 | 
| Rigid microphone arrays                           | ☑    | - 	 | 
| Scattering simulation                             | ☑    | - 	 | 
| Spherical arrays with arbitrary capsule distances | ☑    | - 	 | 
| Custom/measured array IRs into room simulation    | ☑    | - 	 | 
|                                                   |      |  	 | 
|                                                   |      |  	 | 
| masp.spherical_array_processing                   |      |  	 | 
| --------------------------------------------------|      |  	 | 
| A2B theoretical conversion                        | ☑    | - 	 | 
| A2B measurement-based conversion                  | ☑    | - 	 | 
| Beamforming                                       | ☑    | ☑ 	 | 
| Plane-wave decomposition                          | ☑    | - 	 | 
| Nullformer                                        | ☑    | - 	 | 
| Adaptive Beamforming                              | ☑    | - 	 | 
| Adaptive Filtering                                | -    | ☑ 	 | 
| DoA Estimation                                    | ☑    | ☑ 	 | 
| Diffuseness Estimation                            | ☑    | - 	 | 
| Diffuse-field coherence                           | ☑    | - 	 | 
| Blind Source Separation                           | -    | ☑ 	 | 
|                                                   |      |  	 | 
|                                                   |      |  	 | 
| additions?                                        |      |  	 | 
| --------------------------------------------------|      |  	 | 
| source counting?                                  |      |  	 | 