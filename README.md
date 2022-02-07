# Brain-Haemodynamics-Metabolism

The model was implemented in the open source Brain Circulation Model Developer (BCMD) environment, which uses the RADAU5 library to numerically solve the models’ differential-algebraic equation systems. The software and model implementations are available [here](http://www.medphys.ucl.ac.uk/braincirc/). The BCMD environment installation required the gfortran package as well as an updated version of Python (in this report v3.8 was used) and Scipy . More information about the installation requirements and steps can be found [here](https://github.com/multimodalspectroscopy/BCMD/blob/master/doc/manual.pdf). The BCMD environment comes with a GUI (Graphical User Interface) that can be used to run and plot individual simulations. However, due to the nature of the investigation in this report and the high number of simulations needed to run, the GUI was not a feasible option, and instead, several batch scripts written in Python were used to run and analyse the simulations. All the scripts written to run and analyse the simulations get be found [here](https://github.com/miruuna/Brain-Haemodynamics-Metabolism).

The course of a simulation is controlled by an input file with a simple yet not intuitive format, which specifies the simulation time steps and parameter changes as well as the options to reconfigure the the output. The input files are then used together with the compiled model to run batches of simulations. Examples of input files used to run the simulations can be found [here](https://github.com/miruuna/Brain-Haemodynamics-Metabolism).

The written report presenting the methods and the result can be accessed [here](https://github.com/miruuna/Brain-Haemodynamics-Metabolism/blob/main/%5BREPORT%5DBrain_haemodynamics_Serian.pdf).



