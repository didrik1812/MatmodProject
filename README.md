# Modelling reaction-diffusion activity between neurons
This repository contains code for modelling reaction-diffusion activity between neurons and is part of our project work in TMA4195. The code is built upon [Xavier's](https://github.com/xavierr) code, which is available on [bitbucket](https://bitbucket.org/mrst/2022-matmod/src/master/), which again is built upon the [BattMo](https://github.com/BattMoTeam/BattMo) library.

## Installation

1. Download [BattMo](https://github.com/BattMoTeam/BattMo) by cloning their resporitory into your desired folder, i.e run `git clone --recurse-submodules https://github.com/BattMoTeam/BattMo.git`.
2. Clone this repository, i.e `git clone https://github.com/didrik1812/MatmodProject.git`.
3. Add BattMo to MATLAB path by running `startupBattMo` in the MATLAB terminal (while being in the directory where BattMo is cloned to).
4. Add this project to path by running `startup` in the MATLAB terminal (while being in the directory where this repository is cloned to).
5. To check that everything works, try to run the file `runTest3D`.

## About directories and files

The directory `models` contains the model classes which uses BattMo (and MRST) functionality to solve the reaction-diffusion equations. The directory `Examples` contains files which uses the model classes to run the simulation for different cases. The most important files (for experimenting with the simulation) are

* `runTest3D.m`: Runs the simulation on a 3D cylinder
* `runTest2D.m`: Reduces the dimension from 3D to 2D, and runs the simulation on a 2D circle
* `runTest3DGlia.m`: Implements glia-cells in the 3D cylinder.


## Acknowledgements
We want to adress special thanks to [Xavier](https://github.com/xavierr) for very good guidance and debugging with the MATLAB code. 
