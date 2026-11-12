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
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

#include "TRICAL.h"
#include "fit.h"

/*
Language bindings allocate instances as opaque buffers of this size, so make
sure the structure never outgrows it unnoticed.
*/
typedef char _trical_instance_size_check[
    (sizeof(TRICAL_instance_t) <= 2048u) ? 1 : -1];

/*
True if `x` is neither infinite nor NaN. Works on the bit pattern, so that it
cannot be optimised away by compilers told to assume finite arithmetic.
*/
static int _isfinite_f(float x) {
    uint32_t bits;

    memcpy(&bits, &x, sizeof(bits));
    return (bits & 0x7f800000u) != 0x7f800000u;
}

static int _isfinite3_f(const float v[3]) {
    return _isfinite_f(v[0]) && _isfinite_f(v[1]) && _isfinite_f(v[2]);
}

/* Squared distance between two points; may overflow to infinity. */
static float _distance2(const float a[3], const float b[3]) {
    const float dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];

    return dx * dx + dy * dy + dz * dz;
}

/*
Counts `measurement` as one more in the current run of implausible
measurements of the given `kind`. A run is of measurements that could all be
from one sensor: one too far from the first of the run to be that begins a
new run, as does one of another kind. Returns non-zero if the run has by now
come from more than one place -- which a stuck or saturated sensor's does
not.
*/
static int _run_extended(TRICAL_instance_t *instance,
const float measurement[3], unsigned char kind) {
    const float range = TRICAL_RANGE * instance->priv.scale;
    const float distance2 = _distance2(measurement,
                                       instance->priv.reject_anchor);

    if (instance->priv.reject_count == 0 ||
            instance->priv.reject_kind != kind || !_isfinite_f(distance2) ||
            !(distance2 < range * range)) {
        memcpy(instance->priv.reject_anchor, measurement,
               sizeof(instance->priv.reject_anchor));
        instance->priv.reject_count = 0;
        instance->priv.reject_kind = kind;
        instance->priv.reject_moved = 0;
    } else if (!instance->priv.reject_moved) {
        instance->priv.reject_moved =
            distance2 >= instance->priv.gate * instance->priv.gate;
    }
    if (instance->priv.reject_count < 0xffffu) {
        instance->priv.reject_count++;
    }

    return instance->priv.reject_moved;
}

/* True if `measurement` belongs to cell `index`. */
static int _cell_contains(const TRICAL_instance_t *instance,
unsigned int index, const float measurement[3]) {
    float distance2;

    if (instance->priv.cell_count[index] == 0) {
        return 0;
    }

    /* an overflow counts as "no" */
    distance2 = _distance2(measurement, instance->priv.cell_centre[index]);
    return _isfinite_f(distance2) &&
           distance2 < instance->priv.gate * instance->priv.gate;
}

/*
Makes the cell that `measurement` belongs to the current one: one of those
visited recently if it is in one of them, or else a new cell, in place of
the one visited longest ago.
*/
static void _cell_change(TRICAL_instance_t *instance,
const float measurement[3]) {
    float centre[3];
    unsigned short count = 0;
    unsigned int i, found;

    memcpy(centre, measurement, sizeof(centre));
    for (found = 1; found < TRICAL_NUM_CELLS - 1u; found++) {
        if (_cell_contains(instance, found, measurement)) {
            break;
        }
    }
    if (_cell_contains(instance, found, measurement)) {
        memcpy(centre, instance->priv.cell_centre[found], sizeof(centre));
        count = instance->priv.cell_count[found];
    }

    /* the cells before it move back one place, and it goes first */
    for (i = found; i > 0; i--) {
        memcpy(instance->priv.cell_centre[i],
               instance->priv.cell_centre[i - 1u], sizeof(centre));
        instance->priv.cell_count[i] = instance->priv.cell_count[i - 1u];
    }
    memcpy(instance->priv.cell_centre[0], centre, sizeof(centre));
    instance->priv.cell_count[0] = count;
}

/*
Uses `measurement`, which is in the current cell if `in_cell` is non-zero:
adds it to the running sums, unless its cell has had its quota. Returns
non-zero if it was added and is from a different cell from the measurement
before.
*/
static int _measurement_use(TRICAL_instance_t *instance,
const float measurement[3], const float reference_field[3], int in_cell) {
    const float scale = 1.0f / instance->priv.scale;
    float x[3];
    unsigned int i;

    if (!in_cell) {
        _cell_change(instance, measurement);
    }
    if (instance->priv.cell_count[0] >= instance->priv.cell_quota) {
        return 0;
    }
    instance->priv.cell_count[0]++;

    /*
    The cell is centred on the mean of its first few measurements rather
    than on the first one alone -- which, for a sensor at rest, is by
    selection one that noise has carried unusually far.
    */
    if (instance->priv.cell_count[0] <= TRICAL_GATE_AVERAGE) {
        const float weight = 1.0f / (float)instance->priv.cell_count[0];

        for (i = 0; i < 3; i++) {
            instance->priv.cell_centre[0][i] += weight *
                (measurement[i] - instance->priv.cell_centre[0][i]);
        }
    }

    for (i = 0; i < 3; i++) {
        x[i] = (measurement[i] - instance->priv.origin[i]) * scale;
    }
    _trical_moments_accumulate(instance->priv.moments, x);

    /* Reference direction, if there is one. */
    if (reference_field) {
        const float norm2 = reference_field[0] * reference_field[0] +
                            reference_field[1] * reference_field[1] +
                            reference_field[2] * reference_field[2];

        if (norm2 > 0.0f && _isfinite_f(norm2)) {
            const float inverse = 1.0f / sqrtf(norm2);
            trical_wide_t *cross = instance->priv.cross;
            unsigned int j;

            for (i = 0; i < 3; i++) {
                const float f = reference_field[i] * inverse;

                for (j = 0; j < 3; j++) {
                    cross[i * 3 + j] += (trical_wide_t)(f * x[j]);
                }
                cross[9 + i] += (trical_wide_t)f;
            }
            if (instance->priv.reference_count < UINT_MAX) {
                instance->priv.reference_count++;
            }
        }
    }

    if (instance->priv.pending < 0xffffu) {
        instance->priv.pending++;
    }

    return !in_cell;
}

/* True if the caller has written to `state` since TRICAL last did. */
static int _state_written(const TRICAL_instance_t *instance) {
    return memcmp(instance->state, instance->priv.published,
                  sizeof(instance->state)) != 0;
}

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
}

/*
TRICAL_reset:
Discards everything `instance` has learned and resets its calibration
estimate to zero bias and zero scale error. The field norm, measurement noise
and measurement count are not affected.
*/
void TRICAL_reset(TRICAL_instance_t *instance) {
    assert(instance);

    memset(instance->state, 0, sizeof(instance->state));
    memset(&instance->priv, 0, sizeof(instance->priv));
}

/*
TRICAL_norm_set:
Sets the expected field norm (magnitude) of `instance` to `norm`, which must
be between TRICAL_NORM_MIN and TRICAL_NORM_MAX.

What TRICAL learns is the shape of the surface that the raw measurements lie
on, which does not depend on the norm; the norm only fixes the scale of the
calibrated output. So if estimation is under way, a change of norm rescales
(I + D), and leaves the bias alone. `state` is worked out afresh from what
has been learned each time, so it does not matter how often the norm is
changed.
*/
void TRICAL_norm_set(TRICAL_instance_t *instance, float norm) {
    const float previous = instance ? instance->field_norm : 0.0f;

    assert(instance);
    assert(norm >= TRICAL_NORM_MIN && norm <= TRICAL_NORM_MAX);

    /* in a release build, ignore values that make no sense */
    if (!_isfinite_f(norm) || !(norm >= TRICAL_NORM_MIN) ||
            !(norm <= TRICAL_NORM_MAX) || norm == previous) {
        return;
    }

    instance->field_norm = norm;

    /* a norm that the scale estimate cannot follow is no use */
    if (instance->priv.started && !_state_written(instance) &&
            _trical_republish(instance) != 0) {
        instance->field_norm = previous;
    }
}

/*
TRICAL_norm_get:
Returns the expected field norm (magnitude) of `instance`.
*/
float TRICAL_norm_get(const TRICAL_instance_t *instance) {
    assert(instance);

    return instance->field_norm;
}

/*
TRICAL_noise_set:
Sets the standard deviation in measurement supplied to `instance` to `noise`.
*/
void TRICAL_noise_set(TRICAL_instance_t *instance, float noise) {
    assert(instance);
    assert(noise > 0);

    if (!(noise > 0) || !_isfinite_f(noise)) {
        return;
    }

    instance->measurement_noise = noise;
}

/*
TRICAL_noise_get:
Returns the standard deviation in measurements supplied to `instance`.
*/
float TRICAL_noise_get(const TRICAL_instance_t *instance) {
    assert(instance);

    return instance->measurement_noise;
}

/*
TRICAL_noise_estimate_get:
Returns the standard deviation of the measurement noise as estimated by
`instance` from the measurements themselves (strictly: of their distance
from the fitted ellipsoid), in the units of the measurements. Returns zero if
there is no estimate yet.
*/
float TRICAL_noise_estimate_get(const TRICAL_instance_t *instance) {
    assert(instance);

    return instance->priv.noise_estimate * instance->priv.scale;
}

/*
TRICAL_measurement_count_get:
Returns the number of valid measurements previously provided to `instance`
via TRICAL_estimate_update.
*/
unsigned int TRICAL_measurement_count_get(const TRICAL_instance_t *instance) {
    assert(instance);

    return instance->measurement_count;
}

/*
TRICAL_estimate_update
Updates the calibration estimate of `instance` based on the new data in
`measurement`, and (optionally) the expected field direction
`reference_field`. Call this function with each reading you receive from your
sensor.

Per call this does one of three things:
- nothing, if the measurement belongs to the current cell and the cell has
  had its quota (a sensor at rest costs a handful of operations), or if the
  measurement is not plausible;
- adds a measurement or two to the running sums (a few dozen multiply-adds
  each);
- additionally refits the calibration, if the sensor has moved to another
  cell (and the last fit was not just a few measurements ago).

A measurement that leaves the current cell is not used until the next call,
and not at all if the next measurement is back in that cell: the sensor did
not move, the measurement was a glitch.
*/
void TRICAL_estimate_update(TRICAL_instance_t *instance,
const float measurement[3], const float reference_field[3]) {
    float distance2, scale;
    unsigned int i;
    int in_cell, changed = 0;

    assert(instance);
    assert(measurement);

    /*
    If the caller has written a calibration to `state`, restart from it --
    unless it is unusable, in which case put back what was there.
    */
    if (_state_written(instance) && _trical_start(instance, NULL) != 0) {
        memcpy(instance->state, instance->priv.published,
               sizeof(instance->state));
    }

    /*
    Nothing can be done with a measurement that is not a set of numbers; nor
    with any measurement, while the field norm makes no sense (it is a
    public field, so it need not have come through TRICAL_norm_set).
    */
    if (!_isfinite3_f(measurement) || !_isfinite_f(instance->field_norm) ||
            !(instance->field_norm >= TRICAL_NORM_MIN) ||
            !(instance->field_norm <= TRICAL_NORM_MAX)) {
        return;
    }
    if (instance->measurement_count < UINT_MAX) {
        instance->measurement_count++;
    }

    if (!instance->priv.started && _trical_start(instance, NULL) != 0) {
        return;
    }

    /*
    Cells: a measurement further than the gate distance from the centre of
    the current cell is in another cell. One that is not, when the cell has
    had its quota, is of no further interest.
    */
    in_cell = _cell_contains(instance, 0, measurement);
    if (in_cell) {
        /* whatever left the cell has come straight back */
        instance->priv.candidate_state = 0;

        if (instance->priv.cell_count[0] >= instance->priv.cell_quota) {
            instance->priv.reject_count = 0;
            return;
        }
    }

    /*
    Is the measurement plausible? Not if it is too far from the others:
    further from their mean than a few times their spread, in normalised
    coordinates. And not if it is too far from the ellipsoid of the current
    estimate, once there is one that deserves the name.
    */
    scale = 1.0f / instance->priv.scale;
    distance2 = 0.0f;
    for (i = 0; i < 3; i++) {
        const float d = (measurement[i] - instance->priv.origin[i]) * scale -
                        instance->priv.data_mean[i];

        distance2 += d * d;
    }

    if (!_isfinite_f(distance2) ||
            !(distance2 < instance->priv.data_range2)) {
        /*
        Nowhere near anything seen so far: ignore it. A long run of these,
        from more than one place, means that the sensor's characteristics
        have changed completely, so start again where the measurements now
        are, keeping only the scale factors of the current estimate.
        `state` is left alone until there is a new fit to replace it.
        */
        if (_run_extended(instance, measurement, 1) &&
                instance->priv.reject_count >= TRICAL_REJECT_LIMIT &&
                _trical_start(instance, measurement) != 0) {
            instance->priv.reject_count = 0;
        }
        return;
    } else if (instance->priv.fit_gate > 0.0f) {
        const float limit = instance->field_norm * instance->priv.fit_gate;
        float calibrated[3], norm;

        TRICAL_measurement_calibrate(instance, measurement, calibrated);
        norm = sqrtf(calibrated[0] * calibrated[0] +
                     calibrated[1] * calibrated[1] +
                     calibrated[2] * calibrated[2]);

        if (!_isfinite_f(norm) ||
                !(fabsf(norm - instance->field_norm) < limit)) {
            /*
            A few of these are ignored. But if they keep coming, from more
            than one place, it is the estimate that is wrong (the sensor
            has moved somewhere the measurements so far said little about,
            or something about it has changed), and they are used.
            */
            if (!_run_extended(instance, measurement, 2) ||
                    instance->priv.reject_count <= TRICAL_OUTLIER_RUN) {
                return;
            }
        } else {
            instance->priv.reject_count = 0;
        }
    } else {
        instance->priv.reject_count = 0;
    }

    /* The measurement itself may be passed to mean "no reference". */
    if (reference_field == measurement || (reference_field &&
            !_isfinite3_f(reference_field))) {
        reference_field = NULL;
    }

    /*
    A measurement that has left the current cell waits for the next one. If
    one is waiting, this measurement shows that it was no glitch, so use it
    -- which makes its cell the current one.
    */
    if (!in_cell && instance->priv.cell_count[0] != 0) {
        if (instance->priv.candidate_state != 0) {
            changed = _measurement_use(instance, instance->priv.candidate,
                instance->priv.candidate_state == 2 ?
                instance->priv.candidate_reference : NULL, 0);
            instance->priv.candidate_state = 0;

            in_cell = _cell_contains(instance, 0, measurement);
        }
        if (!in_cell) {
            memcpy(instance->priv.candidate, measurement,
                   sizeof(instance->priv.candidate));
            instance->priv.candidate_state = 1;
            if (reference_field) {
                memcpy(instance->priv.candidate_reference, reference_field,
                       sizeof(instance->priv.candidate_reference));
                instance->priv.candidate_state = 2;
            }
        }
    }
    if (in_cell || instance->priv.cell_count[0] == 0) {
        changed |= _measurement_use(instance, measurement, reference_field,
                                    in_cell);
    }

    /*
    Refit when the sensor has moved to another cell, provided there are
    enough measurements to say anything, and enough have been added since
    the last fit to make it worthwhile: to begin with every one is; later, a
    handful.
    */
    if (changed && !(instance->priv.moments[0] <
                     (TRICAL_WIDE_TYPE)TRICAL_REFIT_MIN_COUNT)) {
        const TRICAL_WIDE_TYPE interval = instance->priv.moments[0] /
            (TRICAL_WIDE_TYPE)TRICAL_REFIT_FRACTION;

        if (!((TRICAL_WIDE_TYPE)instance->priv.pending < interval) ||
                instance->priv.pending >= TRICAL_REFIT_INTERVAL) {
            instance->priv.pending = 0;
            _trical_refit(instance);
        }
    }
}

/*
TRICAL_estimate_get
Copies the calibration bias and scale esimates of `instance` to
`bias_estimate` and `scale_estimate` respectively.
*/
void TRICAL_estimate_get(const TRICAL_instance_t *instance,
float bias_estimate[3], float scale_estimate[9]) {
    assert(instance);
    assert(bias_estimate);
    assert(scale_estimate);
    assert(bias_estimate != scale_estimate);

    /* Copy bias estimate from state[0:3] to the destination vector */
    memcpy(bias_estimate, instance->state, 3 * sizeof(float));

    /* Copy the scale estimate to the destination matrix */
    memcpy(scale_estimate, &instance->state[3], 9 * sizeof(float));
}

/*
TRICAL_estimate_get_ext
Same as TRICAL_estimate_get, but additionally copies the bias and scale
estimate variances to `bias_estimate_variance` and `scale_estimate_variance`.
*/
void TRICAL_estimate_get_ext(const TRICAL_instance_t *instance,
float bias_estimate[3], float scale_estimate[9],
float bias_estimate_variance[3], float scale_estimate_variance[9]) {
    float variance[TRICAL_STATE_DIM];

    TRICAL_estimate_get(instance, bias_estimate, scale_estimate);

    /* A bit of paranoia to avoid potential undefined behaviour */
    assert(bias_estimate_variance);
    assert(scale_estimate_variance);
    assert(bias_estimate != bias_estimate_variance);
    assert(bias_estimate != scale_estimate_variance);
    assert(scale_estimate != bias_estimate_variance);
    assert(scale_estimate != scale_estimate_variance);
    assert(bias_estimate_variance != scale_estimate_variance);

    _trical_variance(instance, variance);

    memcpy(bias_estimate_variance, variance, 3 * sizeof(float));
    memcpy(scale_estimate_variance, &variance[3], 9 * sizeof(float));
}

/*
TRICAL_measurement_calibrate
Calibrates `measurement` based on the current calibration estimates, and
copies the result to `calibrated_measurement`. The `measurement` and
`calibrated_measurement` parameters may be pointers to the same vector.

Implements
B' = (I_{3x3} + D)(B - b)

where B' is the calibrated measurement, I_{3x3} is the 3x3 identity matrix,
D is the scale calibration matrix, B is the raw measurement, and b is the bias
vector.

DO NOT pass the calibrated measurement into TRICAL_estimate_update, as it
needs the raw measurement values to work.
*/
void TRICAL_measurement_calibrate(const TRICAL_instance_t *instance,
const float measurement[3], float calibrated_measurement[3]) {
    const float *s;
    float v[3];

    assert(instance);
    assert(measurement);
    assert(calibrated_measurement);

    s = instance->state;

    v[0] = measurement[0] - s[0];
    v[1] = measurement[1] - s[1];
    v[2] = measurement[2] - s[2];

    /* 3x3 matrix multiply */
    calibrated_measurement[0] = v[0] * (s[3] + 1.0f) + v[1] * s[4] +
                                v[2] * s[5];
    calibrated_measurement[1] = v[0] * s[6] + v[1] * (s[7] + 1.0f) +
                                v[2] * s[8];
    calibrated_measurement[2] = v[0] * s[9] + v[1] * s[10] +
                                v[2] * (s[11] + 1.0f);
}
