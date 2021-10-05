# EfficientJetMarcher


## Setup
This project requires the eigen header-only library for linear algebra. Download eigen and set the user environment variable EIGEN_PATH to the eigen-3.x.x folder. Also add this to the Path environment variable if using gcc compiler.

#### Microsoft Visual Studio
The visual studio solution is located in the mvs folder. Make sure EIGEN_PATH is added under EJM "Properties"->"Configuration Properties"->"VC++ Directories"->"Include Directories".

#### GCC
There is a gcc makefile for BASH included in the src directory. If you have the make utility for bash, run the command "make" from within the src directory. The executable main.exe will be created in the parent directory.


