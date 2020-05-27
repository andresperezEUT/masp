Project under development!
__________________

Python port of the following wonderful MATLAB libraries from A. Politis (https://github.com/polarch/)

- shoebox-roomsim
- array-response-simulator
- spherical-array-processing
- spherical-harmonic-transform (partially)



__________________
Comparison with pyroomacoustics


| Feature 	                                        | masp | pra |	
|---	                                            |:---: |:---:|
|                                                   |      |  	 | 
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
| IR estimation?                                    |      |  	 | 
| source counting?                                  |      |  	 | 