
# EfficientJetMarcher
This project is an implementation of the Efficient Jet Marching method for solving for the quasi-potential of 2 dimensional stochastic differential equations,

dX(t)= b(X(t))dt + ε dW(t), 

where W is a 2-dimensional standard Brownian motion ε is a small parameter, and b is a 2D vector field. The method is described in the paper "An efficient jet marcher for computing the quasipotential for 2D SDEs" - by Nicholas Paskal and Maria Cameron -- https://arxiv.org/abs/2109.03424.

## Setup
This project requires the eigen header-only library for linear algebra. Download eigen and set the user environment variable EIGEN_PATH to the eigen-3.x.x folder. Also add this to the Path environment variable if using gcc compiler.

#### Microsoft Visual Studio
The visual studio solution is located in the mvs folder. Make sure EIGEN_PATH is added under EJM "Properties"->"Configuration Properties"->"VC++ Directories"->"Include Directories".

#### GCC
There is a gcc makefile for BASH included in the src directory. If you have the make utility for bash, run the command "make" from within the src directory. The executable main.exe will be created in the parent directory.

## Getting Started
The project startup file is main.cpp. Instructions on how to run EJM to solve for the quasi-potential on a rectangular mesh, as well as a list of the pre-implemented drift functions b(X), are provided at the top of main.cpp.
