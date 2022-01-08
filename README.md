
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

#### Setting the Computational Domain
Currently only rectangular meshes are supported. The rectangular mesh is constructed in  line 78 of *main.cpp* on the rectangle [x<sub>left</sub>, x<sub>right</sub>] x [y<sub>left</sub>, y<sub>right</sub>] with n<sub>x</sub> and n<sub>y</sub> points along each axis, respectively. The implementation assumes that the point attractor lies at the origin, and hence one should have x<sub>left</sub> < 0 < x<sub>right</sub> and y<sub>left</sub> < 0 < y<sub>right</sub>.

#### Setting the Drift Function
Drifts should be set in the function *driftFunctions.cpp*. All supported drifts are implemented as a switch for some master functions. To add a new drift, one should pick an unused integer to use as a case and implement the formulas for that case for the following functions in *driftFunctions.cpp*.
<ol>
  <li> masterDrift -- set drift field <img src="https://render.githubusercontent.com/render/math?math=b:\mathbb{R}^2 \to \mathbb{R}^2"> </li>
  <li> (Optional) masterGradDrift -- set 1st derivative of drift field <img src="https://render.githubusercontent.com/render/math?math=Db:\mathbb{R}^2 \to \mathbb{R}^{2 \times 2}"> </li>
  <li> (Optional) masterHessDrift -- set 2nd derivative of drift field <img src="https://render.githubusercontent.com/render/math?math=D^2b:\mathbb{R}^2 \to \mathbb{R}^{2 \times 2 \times 2}"> </li>
</ol>

If the first and second derivatives are not set, they will default to finite difference calculations. Currently, the coordinate system requires the point attractor to occur at the origin (0,0), so the drifts may need to be adjusted to new coordinates.

The drift functions take a pointer to void argument to contain any additional function parameters. This is cast as a pointer to a particular struct, containing the key for the above function switch, as well as 2 double parameter values (params->a and params->b). If one requires more parameters to specify the drift functions, one should construct a new struct to cast the pointer to, that contains the switch key and any other desired parameters, and adjust all the functions in *driftFunctions.cpp* accordingly.

The struct containing the switch key and the drift parameters is instantiated in *main.cpp* in lines 87-91.

#### (Optional) Setting the Solution Functions

If in addition one has an a-priori solution for the quasi-potential U and its gradient DU -- these can be set in the same manner in the functions **masterSol** and **masterGradSol** in *driftFunctions.cpp*. 

#### Outputs


#### Run Options
Remaining settings are described in detail in file *main.cpp*.

