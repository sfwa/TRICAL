# TRICAL

TRICAL is a UKF-based iterative scale and bias calibration algorithm for
tri-axial field sensors (e.g. magnetometers).

The implementation is on the unscented filter formulation described in
[Real-Time Attitude-Independent Three-Axis Magnetometer Calibration][1];
performance is similar to TWOSTEP but it's less computationally intensive, and
able to provide iterative calibration estimates.

[1]: http://www.acsu.buffalo.edu/~johnc/mag_cal05.pdf


## Usage

TRICAL is configured with an expected field norm (defaulting to 1.0). In the
case of a magnetometer, this would be the magnitude of *B* at its current
location (as output by the WMM, for example).

The input is a sequence of 3-vectors representing the field readings from the
sensor. These are in the same units as the field norm.

The output is a sequence of 3-vectors representing the estimated bias, and
a sequence of 3x3 symmetric matrices representing the estimated scale error.


## Build instructions

Requires `cmake` version 2.8.7 or higher.

Create a build directory outside the source tree, then use cmake to generate
the makefile.

```
mkdir TRICAL_build
cd TRICAL_build
cmake /path/to/TRICAL
```

Now, build the library using the `make` command.


## Testing

The `googletest` library is used for unit testing. To build the unit tests,
use `make unittest`. The unit tests can then be executed by running
`test/unittest` in the build directory.


## Python module installation

Requires `cmake` version 2.8.7 or higher.

Run `python setup.py install` to build the C shared library and install the
Python interface (the `TRICAL` module) in your `site-packages` directory.

Alternatively, just run `pip install https://github.com/sfwa/TRICAL/archive/master.zip#egg=TRICAL-1.0.0`
to download and install.


## Compiling with Texas Instrumets Code Composer Studio 5

Import the root directory of this project (`TRICAL`) into your workspace. CCS
should search all contained files, and find the project files. Complete the
import, and build.
