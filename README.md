# TRICAL

TRICAL is a real-time scale and bias calibration algorithm for tri-axial
field sensors (e.g. magnetometers and accelerometers).

It's designed for embedded use: C99 with no dependencies beyond `libm`, no
dynamic allocation, single-precision arithmetic in the per-sample path, and
less than a kilobyte of state.


## Overview

TRICAL is configured with an expected field norm (defaulting to 1.0). In the
case of a magnetometer, this would be the magnitude of _B_ at its current
location (as output by the WMM, for example).

The input to the calibration process is a sequence of 3-vectors representing
the field readings from the sensor. These are in the same units as the field
norm.

The calibration parameters used are a 3-vector representing the estimated
bias, and a 3x3 matrix representing the estimated scale factor error. A
calibrated reading is `(I + D)(reading - b)`.

These parameters allow a magnetometer to be corrected for hard and soft iron
distortion, sensor bias, sensor scale error, sensor non-orthogonality, and
(if you supply the expected field direction with each reading) misalignment
of the sensor relative to the field. An accelerometer can be corrected for
bias, scale error, non-orthogonality, and misalignment.


## How it works

The readings of a sensor with a bias and scale errors, turned through a
field of constant magnitude, lie on an ellipsoid. TRICAL fits that ellipsoid:

* Each reading used adds to 35 running sums (the moments of the readings up
  to fourth order), which are all that a least-squares fit of a quadric
  surface needs. So the cost per reading is a few dozen multiply-adds; and
  since the sums don't depend on the order of the readings or on the
  estimate at the time, early readings are worth as much as late ones, and
  there's no linearisation to go wrong.
* The fit itself is a 10x10 symmetric eigenproblem, solved whenever the
  sensor has turned far enough to have something new to say (see "cells"
  below). It includes the correction for the bias that measurement noise
  causes in such fits (Taubin's normalisation, with the "hyper-accurate"
  term and renormalisation iteration of Kanatani), and is as accurate as a
  batch maximum-likelihood fit to the same readings.
* Where the readings don't determine the ellipsoid -- the sensor hasn't
  moved much yet, or has only turned about one axis -- a prior holds the
  scale factors, and less firmly the bias, at their initial values. Its
  strength is derived from the measurement noise, which TRICAL estimates
  from the readings; and it is much stronger for as long as the readings
  don't contradict the initial estimate, so that a calibration which is
  already right stays right.
* An ellipsoid determines the calibration only up to a rotation. Without
  further information TRICAL takes the scale matrix to be symmetric. If the
  expected field direction is supplied with the readings, the rotation that
  best aligns the calibrated readings with it is found as well (Wahba's
  problem, in closed form from twelve more running sums).

The details and references are in the comments in `src/fit.c`.

Earlier versions of TRICAL used an unscented Kalman filter for the same
job. The interface is the same, but the behaviour is not; see "Changes from
version 1" below.


## Usage

First, `#include "TRICAL.h"`, and declare a `TRICAL_instance_t` in whatever
scope is appropriate. Call `TRICAL_init(…)` with your TRICAL instance as a
parameter, then (if desired) set the expected field norm and measurement
noise.

The field norm should be the magnitude of the calibrated readings; TRICAL will
scale your measurements to reach that value. The readings themselves should
be in roughly the same units, i.e. have about that magnitude once the bias is
removed. If they aren't (raw ADC counts, say) then tell TRICAL the conversion
by writing an initial scale estimate to the instance's `state` before the
first update -- see "Units" below.

The measurement noise should be something like the standard deviation of your
measurement error. It's not at all critical: TRICAL estimates the noise from
the readings, and only uses your value (limited to between 3% and 5% of the
field norm) until it has done so.

Once you've initialized the instance, you can start giving it sensor readings
to use in estimating the calibration parameters. These readings are passed in
via `TRICAL_estimate_update(…)`, along with the expected field direction if
you have one (pass `NULL` if not); you can get the current calibration
estimate using `TRICAL_estimate_get(…)`.

To apply the current calibration estimate to a measurement, just call
`TRICAL_measurement_calibrate(…)`.

```c
#include "TRICAL.h"

TRICAL_instance_t global_instance;

/* ... */

void your_init_proc(void) {
    /* ... */

    TRICAL_init(&global_instance);
    TRICAL_norm_set(&global_instance, 60.0);
    TRICAL_noise_set(&global_instance, 1.5);
}

void your_sensor_read_proc(void) {
    float sensor_reading[3], expected_field[3];

    /* ... */

    TRICAL_estimate_update(&global_instance, sensor_reading, expected_field);

    /* ... */

    float calibrated_reading[3];
    TRICAL_measurement_calibrate(&global_instance, sensor_reading,
                                 calibrated_reading);

    /* Now use calibrated_reading as an input to your AHRS or whatever */
}
```

### Things worth knowing

**Cells.** TRICAL groups readings into cells about a tenth of the field norm
across (more if the sensor is noisy), uses at most 64 readings from a cell,
and re-fits when the sensor moves to another cell. It remembers the last
four cells, so going back and forth between a few places doesn't keep adding
readings either. So a sensor at rest costs next to nothing and doesn't swamp
the estimate, and a moving one costs one fit per cell, not per reading. It
also means that `state` only changes when the sensor has moved, and can be a
few readings behind.

**Partial coverage.** If the sensor only turns about one axis, the readings
lie on a ring, which fixes some of the calibration parameters and not others.
The calibrated field will have the right magnitude for the orientations the
sensor has been in; elsewhere it depends on how good the initial estimate
was. If that was right (a calibration you saved earlier, say), it stays
right; if not, the estimate improves as soon as the sensor moves somewhere
new.

**The expected field direction** only affects the rotation part of the
calibration. It needn't be accurate reading-by-reading (errors average out),
and only its direction matters, not its magnitude.

**Saving and restoring.** `state` (12 floats: the bias, then the scale
matrix in row-major order) is the calibration. You can copy it out, and you
can write a saved calibration back into it at any time, including before the
first update: TRICAL notices, and carries on from the calibration you wrote,
treating it as its initial estimate. It is left exactly as you wrote it
until the sensor has moved enough for there to be something to learn, and
it is only changed much if the readings then contradict it. Values that
can't be a calibration (not finite, or a singular scale matrix) are
discarded.

**Units.** TRICAL works out the scale of the readings for itself if they are
between about a quarter and six times the size that the field norm implies
(more slowly, the further off they are). If they are in quite different
units from the field norm, say so by writing the conversion as the initial
scale estimate, after `TRICAL_init(…)` and before the first update. For a
sensor that reports about 1500 counts per unit of field:

```c
global_instance.state[3] = global_instance.state[7] =
    global_instance.state[11] = 1.0f / 1500.0f - 1.0f;
```

(`state` holds the scale matrix less the identity in single precision, so a
conversion factor of 1/1500 is held to about one part in 10,000, and a
smaller one to less. If that matters, change the units of the field norm
instead.)

**The field norm** can be changed at any time, and as often as you like --
if it comes from a field model that is evaluated for every reading, say.

**Starting again.** `TRICAL_reset(…)` discards what has been learned and
returns the estimate to zero bias and no scale error (so write your initial
scale estimate again, if you used one), keeping the field norm and noise
settings. Use it if the sensor or its surroundings change: TRICAL doesn't
forget old readings by itself. (It does start again by itself after a few
hundred readings in a row that are nowhere near the old ones, if they hang
together and aren't all the same.)

**How good is it?** `TRICAL_noise_estimate_get(…)` returns the RMS distance
of the readings from the fitted ellipsoid, in the units of the readings. If
that's much more than your sensor's noise, the readings don't lie on an
ellipsoid -- there's interference, or the field norm is changing.
`TRICAL_estimate_get_ext(…)` gives approximate variances of the estimates,
which are large for whatever the readings so far don't determine.

**Bad readings.** These are not used:

* readings that aren't finite;
* readings several field norms away from all the others;
* once there is an estimate to judge by, readings whose calibrated magnitude
  is more than 25% (or 8 times the noise, if that is more) from the field
  norm -- unless 64 of them arrive in a row from more than one place, which
  means that it's the estimate that is wrong;
* a reading that is somewhere else entirely when the next one is back where
  the sensor was: that's a glitch, not a movement. (To tell, a reading that
  shows the sensor to have moved isn't used until the next one arrives.)

So a sensor that is stuck, saturated or disconnected for a while does no
harm, however long it lasts.


## Resource usage

Measured for version 2.0 with gcc; the Cortex-M4F figures are for
`-mcpu=cortex-m4 -mfpu=fpv4-sp-d16 -mfloat-abi=hard -Os`, and the stack
figures are sums of the frame sizes the compiler reports along the deepest
call chain.

| | default | `TRICAL_WIDE_TYPE=float` |
|---|---|---|
| Instance size | 840 bytes | 568 bytes |
| Code size (Cortex-M4F) | 19 kB | 14 kB |
| Stack in `TRICAL_estimate_update` (Cortex-M4F) | 3.1 kB | 1.6 kB |
| Stack in `TRICAL_estimate_get_ext` (Cortex-M4F) | 3.4 kB | 2.0 kB |
| Stack in `TRICAL_norm_set` (Cortex-M4F) | 1.2 kB | 0.7 kB |

`TRICAL_estimate_update` does one of three things, with very different costs
(x86-64 instructions, and approximate time on a 2.1 GHz Xeon):

| | instructions | time |
|---|---|---|
| Sensor at rest: nothing to do | 120 | 0.01 µs |
| Reading added to the running sums | 400-700 (550-850 with a reference) | 0.05 µs |
| Reading added, and the estimate refitted | 21,000 (23,000 with a reference) | 1.9 µs |
| Average, sensor turning 1° per reading | 2,400 | 0.25 µs |
| Average, sensor turning 10° per reading | 3,300 | 0.3 µs |

A refit happens when the sensor has turned into a new cell, and at most once
every 8 readings (once there are more than a hundred or so). It normally
needs one matrix factorisation; readings that are hard to make sense of can
take more, up to a fixed limit of 24.

`TRICAL_measurement_calibrate` takes about 36 instructions.

The running sums and the fit use double precision by default. Everything that
happens for every reading is single precision apart from the additions to
the sums (35, or 47 with a reference), so that's a modest cost even without a
double-precision FPU; the refit is where double precision costs most. If you
need to avoid double precision altogether, build with
`-DTRICAL_WIDE_TYPE=float` (everything that includes `TRICAL.h` must be built
the same way). The calibration is then accurate to about 0.01% of the field
norm at best, rather than limited only by the noise of the sensor; each time
131,072 readings have been used the weight of them all is halved, so the
oldest gradually count for less; and where the readings cover only a small
part of the ellipsoid, what they leave open is held less precisely (to a
percent or two of the field norm).

If the time or stack that a refit takes is a problem where your sensor is
read, call `TRICAL_estimate_update` from somewhere less critical: nothing
depends on its being called at regular intervals, or for every reading.


## Build instructions

Requires `cmake` version 3.10 or higher (or just compile `src/TRICAL.c` and
`src/fit.c` with `include` and `src` on the include path).

Create a build directory outside the source tree, then use cmake to generate
the makefile.

```
mkdir TRICAL_build
cd TRICAL_build
cmake /path/to/TRICAL
```

Now, build the library using the `make` command.


## Testing

The `googletest` library is used for unit testing; an installed copy is used
if there is one, otherwise cmake (3.14 or higher) downloads it. To build the
unit tests, use `make unittest`. The unit tests can then be executed by
running `test/unittest` in the build directory.

`make benchmark` builds `test/benchmark`, which times the library on
simulated data.


## Python module installation

Requires `cmake` version 3.10 or higher, and Python 3.

Run `pip install .` in this directory to build the C shared library and
install the Python interface (the `TRICAL` module).

Alternatively, build the shared library as above and set the `TRICAL_LIBRARY`
environment variable to its path; the module can then be used straight from
the `python` directory.

Once the module has been installed, you can run the TRICAL batch calibration
script as follows:

```
python -m TRICAL 1.0 1e-6 < /path/to/input.csv > /path/to/output.csv
```

(Where `1.0` is the field norm, and `1e-6` is the measurement noise. Change
those values as appropriate.)

The script reads 3 comma-separated input values per `\n`-terminated line from
`stdin`, and writes 3 comma-separated calibrated output values to `stdout` per
input line read. If a line has 6 values, the last 3 are taken to be the
expected field. Upon completion (EOF on `stdin`), the script outputs the
calibration estimate to `stderr`.

You can also generate a WebGL point cloud visualisation of the data by
including `html` after the command:

```
python -m TRICAL 1.0 1e-6 html < /path/to/input.csv > /path/to/output.html
```

The HTML data is written to `stdout`, and can be opened in any
WebGL-supporting browser.


## Compiling with Texas Instruments Code Composer Studio 5

Import the root directory of this project (`TRICAL`) into your workspace. CCS
should search all contained files, and find the project files. Complete the
import, and build.


## Changes from version 1

The functions and the layout of `state` are as they were, so existing code
should compile and run unchanged. What's different:

* `TRICAL_estimate_update(…)` no longer requires the expected field; pass
  `NULL` (or the reading itself, as some callers did) if you don't have it.
  With the expected field, the calibration includes the rotation of the
  sensor relative to it; without, the scale matrix is symmetric.
* The estimate changes when the sensor has moved to another cell, rather
  than with every reading; and readings that don't look right are not used
  (see "Bad readings" above).
* The measurement noise setting is no longer a tuning parameter, and the
  field norm no longer needs to be close to 1 for good convergence.
* `TRICAL_instance_t` has different private contents (and no
  `state_covariance`); it's 840 rather than 636 bytes.
* An update is much cheaper on average (see above), but an update that
  refits the estimate takes about three times as long as every update used
  to, and needs about 3 kB of stack rather than 1 kB. The library is about
  19 kB of code on a Cortex-M4F rather than 3 kB.
* New functions: `TRICAL_reset(…)` and `TRICAL_noise_estimate_get(…)`.
* `state` can be written as well as read.
