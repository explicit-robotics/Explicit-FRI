# Expplicit-FRI
This repository contains Exp[licit]-FRI -- a software to control robots via the KUKA FRI interface (C++). 

The current version of this software is in active development, i.e., might be unstable. A formal announcement will be made once the first stable version of Exp[licit]-cpp is complete.

# Library structure
While any folder structure can be created by adapting the CMAKE-files, this is the initial setting:

Explicit-FRI (head folder)<br />
&nbsp;&nbsp;&nbsp;|--- Libraries<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|--- Explicit-cpp<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|--- Eigen<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|--- myFRIClient<br />
&nbsp;&nbsp;&nbsp;|--- ClientApplications<br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|--- basic_app

-*Eigen*: To run the code, you will have to install [Eigen](https://eigen.tuxfamily.org/dox/GettingStarted.html) on your computer and link the directory (usually /usr/include/) in the CMAKE files of *Explicit-cpp* and *myFRIClient*

-*Explicit-cpp*: All kinematics and dynamics are calculated with Exp[licit]-cpp, a robotic modelling software based on differential geometry. The source files and can be [found here](https://github.com/explicit-robotics/Explicit-cpp). This software contains the submodule Exp[licit]-cpp. Exp[licit]-cpp in turn contains [Eigen version 3.4.0](https://gitlab.com/libeigen/eigen/-/releases/3.4.0). To update the submodules, type:

```
    git submodule update --init --recursive
```
We always need to check whether Explicit-cpp is updated. For that, 
```
    git submodule update --remote
```

To run Exp[licit]-FRI, a .so-file has to be compiled. This can be done by running

```
    make -f Makefile-lib
```

 Note that the source files are created for a "Media Flange Touch." If needed, please adapt the Flange Position in *exp_robots.cpp*.  

-*myFRIClient*: This library is provided by KUKA and establishes a state-machine for communication between the Client-PC and the Sunrise controller. The software uploaded here was extracted from a Sunrise 1.17-version. For compilation, please refer to the README-file inside the folder.

-*basic_app*: To get started with FRI (C++), we provide a basic application. The application was compiled by using [Qt-Creator](https://www.qt.io/product/development-tools). Simply run the CMAKE-file to compile the application.

The source file includes an iir-filter, with coefficients determined [here](http://www.winfilter.20m.com/). The filter is needed to activate the robot's build-in friction compensation. Before sending the torques, we add a simple mean filter to smooth out the torque signals. The basic application calculates the Forward Kinematics, Jacobian Matrix, and Mass matrix of the robot and prints the calculational effort. 

# Authors
This software was developed by [Moses Chong-ook Nah](https://mosesnah-shared.github.io/) and [Johannes Lachner](https://jlachner.github.io/), working in the Newman Laboratory for Biomechanics and Human Rehabilitation at MIT.

# Documentation of the Software 
Detailed documentation of the software can be [checked here](https://explicit-robotics.github.io/).

# Expplicit-MATLAB
If you are interested in simulating the robot first, check out [Explicit-MATLAB](https://github.com/explicit-robotics/Explicit-MATLAB).
