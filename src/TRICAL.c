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

#include <assert.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include "TRICAL.h"
#include "filter.h"

/*
TRICAL_init:
Initializes `instance`. Must be called prior to any other TRICAL procedures
taking `instance` as a parameter.

If called on an instance that has already been initialized, TRICAL_init will
reset that instance to its default state.
*/
void TRICAL_init(TRICAL_instance_t *instance) {
    assert(instance);

    memset(instance, 0, sizeof(TRICAL_instance_t));

    instance->field_norm = 1.0f;
    instance->measurement_noise = 1e-6f;

    /*
    Set the state covariance diagonal to a small value, so that we can run the
    Cholesky decomposition without blowing up
    */
    unsigned int i;
    for (i = 0; i < TRICAL_STATE_DIM * TRICAL_STATE_DIM; i += 10) {
        instance->state_covariance[i] = 1e-2f;
    }
}

/*
TRICAL_norm_set:
Sets the expected field norm (magnitude) of `instance` to `norm`. If `norm`
differs from the instance's current field norm, all estimates are multiplied
by the ratio of new norm:old norm.
*/
void TRICAL_norm_set(TRICAL_instance_t *instance, float norm) {
    assert(instance);
    assert(norm > 0.0f);

    if (fabs(norm - instance->field_norm) < FLT_EPSILON) {
        return;
    }

    /*
    Since the new norm differs from the current norm, we need to re-scale our
    state and covariance estimates accordingly.
    */
    float scale = norm / instance->field_norm;
    unsigned int i;

    #pragma MUST_ITERATE(9, 9)
    for (i = 0; i < TRICAL_STATE_DIM; i++) {
        instance->state[i] *= scale;
    }

    scale = (norm * norm) / (instance->field_norm * instance->field_norm);

    #pragma MUST_ITERATE(81, 81)
    for (i = 0; i < TRICAL_STATE_DIM * TRICAL_STATE_DIM; i++) {
        instance->state_covariance[i] *= scale;
    }

    instance->field_norm = norm;
}

/*
TRICAL_norm_get:
Returns the expected field norm (magnitude) of `instance`.
*/
float TRICAL_norm_get(TRICAL_instance_t *instance) {
    assert(instance);

    return instance->field_norm;
}

/*
TRICAL_noise_set:
Sets the standard deviation in measurement supplied to `instance` to `noise`.
*/
void TRICAL_noise_set(TRICAL_instance_t *instance, float noise) {
    assert(instance);
    assert(noise > 0.0f);

    instance->measurement_noise = noise;
}

/*
TRICAL_noise_get:
Returns the standard deviation in measurements supplied to `instance`.
*/
float TRICAL_noise_get(TRICAL_instance_t *instance) {
    assert(instance);

    return instance->measurement_noise;
}

/*
TRICAL_measurement_count_get:
Returns the number of measurements previously provided to `instance` via
TRICAL_estimate_update.
*/
unsigned int TRICAL_measurement_count_get(TRICAL_instance_t *instance) {
    assert(instance);

    return instance->measurement_count;
}

/*
TRICAL_estimate_update
Updates the calibration estimate of `instance` based on the new data in
`measurement`. Call this function with each reading you receive from your
sensor.
*/
void TRICAL_estimate_update(TRICAL_instance_t *instance,
float measurement[3]) {
    assert(instance);
    assert(measurement);

    _trical_filter_iterate(instance, measurement);
    instance->measurement_count++;
}

/*
TRICAL_estimate_get
Copies the calibration bias and scale esimates of `instance` to
`bias_estimate` and `scale_estimate` respectively. A new calibration estimate
will be available after every call to TRICAL_estimate_update.
*/
void TRICAL_estimate_get(TRICAL_instance_t *restrict instance,
float bias_estimate[3], float scale_estimate[9]) {
    assert(instance);
    assert(bias_estimate);
    assert(scale_estimate);
    assert(bias_estimate != scale_estimate);

    /* Copy bias estimate from state[0:3] to the destination vector */
    memcpy(bias_estimate, instance->state, 3 * sizeof(float));

    /*
    Now copy scale estimate/scale estimate covariance. Since the scale matrix
    is symmetric, we only store the upper triangle values in the state vector,
    like this:
    3  4  5
       6  7
          8

    So, on output, we need to copy values 4, 5 and 7 to two locations:

    3  4  5
    4  6  7
    5  7  8
    */
    scale_estimate[0] = instance->state[3];
    scale_estimate[1] = instance->state[4];
    scale_estimate[2] = instance->state[5];

    scale_estimate[3] = instance->state[4];
    scale_estimate[4] = instance->state[6];
    scale_estimate[5] = instance->state[7];

    scale_estimate[6] = instance->state[5];
    scale_estimate[7] = instance->state[7];
    scale_estimate[8] = instance->state[8];
}

/*
TRICAL_estimate_get_ext
Same as TRICAL_estimate_get, but additionally copies the bias and scale
estimate variances to `bias_estimate_variance` and `scale_estimate_variance`.
*/
void TRICAL_estimate_get_ext(TRICAL_instance_t *restrict instance,
float bias_estimate[3], float scale_estimate[9],
float bias_estimate_variance[3], float scale_estimate_variance[9]) {
    TRICAL_estimate_get(instance, bias_estimate, scale_estimate);

    /* A bit of paranoia to avoid potential undefined behaviour */
    assert(bias_estimate_variance);
    assert(scale_estimate_variance);
    assert(bias_estimate != bias_estimate_variance);
    assert(bias_estimate != scale_estimate_variance);
    assert(scale_estimate != bias_estimate_variance);
    assert(scale_estimate != scale_estimate_variance);
    assert(bias_estimate_variance != scale_estimate_variance);

    /* Copy bias estimate covariance from the state covariance diagonal */
    bias_estimate_variance[0] = instance->state_covariance[0 * 9 + 0];
    bias_estimate_variance[1] = instance->state_covariance[1 * 9 + 1];
    bias_estimate_variance[2] = instance->state_covariance[2 * 9 + 2];

    /*
    Now copy scale estimate covariance. This is the same approach as we used
    to extract the scale estimate itself (see TRICAL_estimate_get), except
    we're pulling the values from the diagonal of the state covariance matrix
    so the indices are:

    3*9+3    4*9+4    5*9+5
    4*9+4    6*9+6    7*9+7
    5*9+5    7*9+7    8*9+8
    */
    scale_estimate_variance[0] = instance->state_covariance[3 * 9 + 3];
    scale_estimate_variance[1] = instance->state_covariance[4 * 9 + 4];
    scale_estimate_variance[2] = instance->state_covariance[5 * 9 + 5];

    scale_estimate_variance[3] = instance->state_covariance[4 * 9 + 4];
    scale_estimate_variance[4] = instance->state_covariance[6 * 9 + 6];
    scale_estimate_variance[5] = instance->state_covariance[7 * 9 + 7];

    scale_estimate_variance[6] = instance->state_covariance[5 * 9 + 5];
    scale_estimate_variance[7] = instance->state_covariance[7 * 9 + 7];
    scale_estimate_variance[8] = instance->state_covariance[8 * 9 + 8];
}

/*
TRICAL_measurement_calibrate
Calibrates `measurement` based on the current calibration estimates, and
copies the result to `calibrated_measurement`. The `measurement` and
`calibrated_measurement` parameters may be pointers to the same vector.

DO NOT pass the calibrated measurement into TRICAL_estimate_update, as it
needs the raw measurement values to work.
*/
void TRICAL_measurement_calibrate(TRICAL_instance_t *restrict instance,
float measurement[3], float calibrated_measurement[3]) {
    assert(instance);
    assert(measurement);
    assert(calibrated_measurement);

    /* Pass off to the internal function */
    _trical_measurement_calibrate(instance->state, measurement,
                                  calibrated_measurement);
}
