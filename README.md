# Expplicit-FRI
This repository contains Exp[licit]-FRI -- a software to control robots via the KUKA FRI interface (C++). 

This software contains the submodule Exp[licit]-cpp. Exp[licit]-cpp in turn contains [Eigen version 3.4.0](https://gitlab.com/libeigen/eigen/-/releases/3.4.0). To update the submodules, type:
```
    git submodule update --init --recursive
```

We always need to check whether Explicit-cpp is updated. For that, 
```
    git submodule update --remote
```

The current version of this software is in active development, i.e., might be unstable. A formal announcement will be made once the first stable version of Exp[licit]-cpp is complete.

# Authors
This software was developed by [Moses Chong-ook Nah](https://mosesnah-shared.github.io/) and [Johannes Lachner](https://jlachner.github.io/), working in the Newman Laboratory for Biomechanics and Human Rehabilitation at MIT.

# Documentation of the Software 
Detailed documentation of the software can be [checked here](https://explicit-robotics.github.io/).

# Expplicit-MATLAB
If you are interested in simulating the robot first, check out [Explicit-MATLAB](https://github.com/explicit-robotics/Explicit-MATLAB).
