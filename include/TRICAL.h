/*
Copyright (C) 2013 Ben Dyer

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef TRICAL_H
#define TRICAL_H

#ifdef __cplusplus
extern "C" {
#endif

/*
Number of calibration parameters: a 3-vector of biases followed by a 3x3
scale matrix in row-major order.
*/
#define TRICAL_STATE_DIM 12

/* The range of field norms that TRICAL accepts. */
#define TRICAL_NORM_MIN (1e-12f)
#define TRICAL_NORM_MAX (1e12f)

/*
Number of running sums kept by the estimator: the sums of all monomials
x^a y^b z^c of the (normalised) measurements with a + b + c <= 4.
*/
#define TRICAL_NUM_MOMENTS 35

/*
Number of cells of measurements that TRICAL remembers having visited: a
sensor that goes back and forth between this many places or fewer does not
keep adding to the weight of the measurements made there.
*/
#define TRICAL_NUM_CELLS 4

/*
Floating-point type used for the running sums and for the fit. The sensor
hot path (TRICAL_measurement_calibrate, and the per-sample part of
TRICAL_estimate_update) is always single precision.

`double` is strongly recommended: the fit extracts a quantity of the order
of the squared noise-to-signal ratio from sums of order one, which needs
more than the 24 bits of a `float` unless the sensor is noisy and the bias
small. Defining TRICAL_WIDE_TYPE as `float` is supported for targets where
double-precision arithmetic is unaffordable, at reduced accuracy. The
library and all code that uses TRICAL_instance_t must agree on it.
*/
#ifndef TRICAL_WIDE_TYPE
#define TRICAL_WIDE_TYPE double
#endif

typedef struct {
    /* Expected magnitude of a calibrated measurement. */
    float field_norm;

    /* Measurement noise assumed until it can be estimated from the data. */
    float measurement_noise;

    /*
    The calibration estimate: the bias b in state[0..2], and the scale matrix
    D in state[3..11] (row-major), such that

        calibrated = (I + D) (measurement - b)

    This array is part of the public interface. It may be read at any time,
    and it may be written, e.g. to restore a calibration saved earlier:
    estimation then restarts from the values written, which are used as the
    initial estimate.

    Note that D is the difference between the scale matrix and the identity,
    in single precision. A scale matrix much smaller than the identity (raw
    measurements in units much smaller than those of the field norm) is
    therefore held with correspondingly less relative precision: about 1e-7
    divided by its size.
    */
    float state[TRICAL_STATE_DIM];

    /* Number of valid measurements passed to TRICAL_estimate_update. */
    unsigned int measurement_count;

    /*
    Everything below is private to the implementation. It is declared here
    only so that instances can be allocated statically.
    */
    struct {
        /* running sums of monomials of the normalised measurements */
        TRICAL_WIDE_TYPE moments[TRICAL_NUM_MOMENTS];
        /* running sums of (reference direction) x (normalised measurement) */
        TRICAL_WIDE_TYPE cross[12];
        /* current and prior quadric (ellipsoid) coefficients */
        TRICAL_WIDE_TYPE quadric[10];
        TRICAL_WIDE_TYPE quadric_prior[10];

        /* `state` as last written by TRICAL: detects writes by the caller */
        float published[TRICAL_STATE_DIM];
        /* rotation part of the initial estimate */
        float rotation_prior[9];
        /* origin and scale of the normalised measurement coordinates */
        float origin[3];
        float scale;
        /*
        centres of the cells of measurements visited most recently (the
        current one first), and the size of a cell
        */
        float cell_centre[TRICAL_NUM_CELLS][3];
        float gate;
        /*
        a measurement that left the current cell, and its reference field,
        held until the next one shows whether the sensor really moved
        */
        float candidate[3];
        float candidate_reference[3];
        /* mean of the data, and squared range around it that is accepted */
        float data_mean[3];
        float data_range2;
        /*
        estimated measurement noise in normalised units; zero until an
        estimate is available
        */
        float noise_estimate;
        /* diagonal shift that the last fit needed */
        float shift;

        /*
        how far, as a fraction of the field norm, the magnitude of a
        calibrated measurement may be from it; zero until there is an
        estimate to judge by
        */
        float fit_gate;
        /* where the current run of implausible measurements began */
        float reject_anchor[3];

        unsigned int reference_count;
        /*
        measurements used from each of those cells (zero: no such cell), and
        the most allowed
        */
        unsigned short cell_count[TRICAL_NUM_CELLS];
        unsigned short cell_quota;
        /* measurements used since the last fit */
        unsigned short pending;
        /* length of the current run of implausible measurements */
        unsigned short reject_count;
        unsigned char started;
        /* 0: no candidate; 1: candidate; 2: candidate with reference */
        unsigned char candidate_state;
        /*
        the current run: what makes its measurements implausible (1: nowhere
        near the others; 2: not on the ellipsoid of the estimate), and
        whether it has come from more than one place
        */
        unsigned char reject_kind;
        unsigned char reject_moved;
    } priv;
} TRICAL_instance_t;

/*
TRICAL_init:
Initializes `instance`. Must be called prior to any other TRICAL procedures
taking `instance` as a parameter.

If called on an instance that has already been initialized, TRICAL_init will
reset that instance to its default state: a field norm of 1.0, a measurement
noise of 1e-6, zero bias and zero scale error.
*/
void TRICAL_init(TRICAL_instance_t *instance);

/*
TRICAL_reset:
Discards everything `instance` has learned and resets its calibration
estimate to zero bias and zero scale error. The field norm, measurement noise
and measurement count are not affected.
*/
void TRICAL_reset(TRICAL_instance_t *instance);

/*
TRICAL_norm_set:
Sets the expected field norm (magnitude) of `instance` to `norm`, which must
be between TRICAL_NORM_MIN and TRICAL_NORM_MAX. Calibrated measurements are
scaled to have this magnitude.

Set the norm before supplying measurements, and supply measurements in the
same units as the norm. If the norm is changed once estimation is under way,
the scale estimate follows it (it is worked out afresh from what has been
learned, so the norm can be changed as often as you like); the bias estimate
is not affected.

A norm outside the range, or one for which the scale estimate could not be
represented, is ignored.
*/
void TRICAL_norm_set(TRICAL_instance_t *instance, float norm);

/*
TRICAL_norm_get:
Returns the expected field norm (magnitude) of `instance`.
*/
float TRICAL_norm_get(const TRICAL_instance_t *instance);

/*
TRICAL_noise_set:
Sets the standard deviation of the noise in measurements supplied to
`instance` to `noise` (in the same units as the measurements). This must be
positive.

This is only an initial assumption: once enough measurements have been
received, TRICAL estimates the noise level from the data and uses that
instead (see TRICAL_noise_estimate_get). Until then it assumes the value set
here, limited to between 3% and 5% of the field magnitude. The assumed noise
level determines how strongly the estimate is held to its initial value where
the measurements received so far do not determine it.
*/
void TRICAL_noise_set(TRICAL_instance_t *instance, float noise);

/*
TRICAL_noise_get:
Returns the value set by TRICAL_noise_set.
*/
float TRICAL_noise_get(const TRICAL_instance_t *instance);

/*
TRICAL_noise_estimate_get:
Returns the standard deviation of the measurement noise as estimated by
`instance` from the measurements themselves (strictly: of their distance
from the fitted ellipsoid), in the units of the measurements. Returns zero if
there is no estimate yet.

This is a useful check on the quality of a calibration: a value much larger
than the noise of the sensor means that the measurements do not lie on an
ellipsoid, e.g. because of interference that varies with time.
*/
float TRICAL_noise_estimate_get(const TRICAL_instance_t *instance);

/*
TRICAL_measurement_count_get:
Returns the number of valid (finite) measurements previously provided to
`instance` via TRICAL_estimate_update.
*/
unsigned int TRICAL_measurement_count_get(const TRICAL_instance_t *instance);

/*
TRICAL_estimate_update
Updates the calibration estimate of `instance` based on the new data in
`measurement`, the raw sensor reading. Call this function with each reading
you receive from your sensor.

`reference_field` is optional. It is the expected direction of the field in
the frame in which calibrated measurements are wanted (e.g. the body frame),
at the time of the measurement; only its direction is used. If it is given,
the scale estimate additionally includes the rotation that best aligns
calibrated measurements with the reference. Pass NULL (or `measurement`
itself) if no reference is available: the scale estimate is then symmetric,
i.e. free of rotation.

The bias and the symmetric part of the scale estimate are determined from the
measurements alone, and are not affected by errors in the reference.

Measurements that are implausible are not used: ones that are not finite;
ones that are nowhere near the others, or not near the ellipsoid of the
current estimate (unless they keep coming, in which case it is the estimate
that needs to change); and one that is somewhere else when the next is back
where the sensor was before. To tell the last kind, a measurement that shows
the sensor to have moved is not used until the next call.

The call does a bounded amount of work: usually a few dozen arithmetic
operations (next to none for a sensor at rest), occasionally -- when the
sensor has turned far enough for its measurements to have something new to
say -- a refit of the estimate, taking a few thousand.
*/
void TRICAL_estimate_update(TRICAL_instance_t *instance,
const float measurement[3], const float reference_field[3]);

/*
TRICAL_estimate_get
Copies the calibration bias and scale estimates of `instance` to
`bias_estimate` and `scale_estimate` respectively.
*/
void TRICAL_estimate_get(const TRICAL_instance_t *instance,
float bias_estimate[3], float scale_estimate[9]);

/*
TRICAL_estimate_get_ext
Same as TRICAL_estimate_get, but additionally copies approximate variances of
the bias and scale estimates to `bias_estimate_variance` and
`scale_estimate_variance`. The variance of the scale estimate does not
include the uncertainty of the rotation derived from reference directions.

This is considerably more expensive than TRICAL_estimate_get.
*/
void TRICAL_estimate_get_ext(const TRICAL_instance_t *instance,
float bias_estimate[3], float scale_estimate[9],
float bias_estimate_variance[3], float scale_estimate_variance[9]);

/*
TRICAL_measurement_calibrate
Calibrates `measurement` based on the current calibration estimates, and
copies the result to `calibrated_measurement`. The two parameters may point
to the same vector.

DO NOT pass the calibrated measurement into TRICAL_estimate_update, as it
needs the raw measurement values to work.
*/
void TRICAL_measurement_calibrate(const TRICAL_instance_t *instance,
const float measurement[3], float calibrated_measurement[3]);

#ifdef __cplusplus
}
#endif

#endif
