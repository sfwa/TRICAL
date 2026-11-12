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

#include <gtest/gtest.h>

#include <cmath>
#include <cstring>
#include <limits>

#include "TRICAL.h"
#include "support.h"

using support::Random;
using support::Sensor;
using support::Tumble;
using support::Vec3;

namespace {

/* Feeds `count` readings of a smoothly tumbling `sensor` to `cal`. */
void tumble(TRICAL_instance_t &cal, const Sensor &sensor, Random &rng,
        int count, bool with_reference = false) {
    Tumble motion(rng);

    for (int i = 0; i < count; i++) {
        const Vec3 direction = motion.step(0.03);
        const Vec3 reading = sensor.read(direction, rng);

        support::update(cal, reading, with_reference ? &direction : nullptr);
    }
}

} /* namespace */

TEST(TRICAL, Initialisation) {
    TRICAL_instance_t cal;

    /* make sure init does not depend on what was there before */
    std::memset(&cal, 0xa5, sizeof(cal));

    TRICAL_init(&cal);
    EXPECT_FLOAT_EQ(1.0f, TRICAL_norm_get(&cal));
    EXPECT_FLOAT_EQ(1e-6f, TRICAL_noise_get(&cal));
    EXPECT_EQ(0u, TRICAL_measurement_count_get(&cal));

    float bias[3], scale[9];
    TRICAL_estimate_get(&cal, bias, scale);
    for (int i = 0; i < 3; i++) {
        EXPECT_EQ(0.0f, bias[i]);
    }
    for (int i = 0; i < 9; i++) {
        EXPECT_EQ(0.0f, scale[i]);
    }
}

/* Check that the field norm can be set and read */
TEST(TRICAL, NormGetSet) {
    TRICAL_instance_t cal;

    TRICAL_init(&cal);
    TRICAL_norm_set(&cal, 2.0f);
    EXPECT_FLOAT_EQ(2.0f, TRICAL_norm_get(&cal));
}

/* Check that the measurement noise can be set and read */
TEST(TRICAL, NoiseGetSet) {
    TRICAL_instance_t cal;

    TRICAL_init(&cal);
    TRICAL_noise_set(&cal, 2.0f);
    EXPECT_FLOAT_EQ(2.0f, TRICAL_noise_get(&cal));
}

#ifdef NDEBUG
/* With assertions off, nonsensical settings are ignored */
TEST(TRICAL, InvalidSettingsAreIgnored) {
    const float nan = std::numeric_limits<float>::quiet_NaN();
    const float inf = std::numeric_limits<float>::infinity();
    TRICAL_instance_t cal;

    TRICAL_init(&cal);
    TRICAL_norm_set(&cal, 3.0f);
    TRICAL_noise_set(&cal, 0.5f);

    TRICAL_norm_set(&cal, 0.0f);
    TRICAL_norm_set(&cal, -1.0f);
    TRICAL_norm_set(&cal, nan);
    TRICAL_norm_set(&cal, inf);
    TRICAL_norm_set(&cal, TRICAL_NORM_MIN / 10.0f);
    TRICAL_norm_set(&cal, TRICAL_NORM_MAX * 10.0f);
    EXPECT_FLOAT_EQ(3.0f, TRICAL_norm_get(&cal));

    TRICAL_noise_set(&cal, 0.0f);
    TRICAL_noise_set(&cal, -1.0f);
    TRICAL_noise_set(&cal, nan);
    TRICAL_noise_set(&cal, inf);
    EXPECT_FLOAT_EQ(0.5f, TRICAL_noise_get(&cal));
}
#endif

/* Check that the measurement count can be read */
TEST(TRICAL, MeasurementCountGet) {
    TRICAL_instance_t cal;

    TRICAL_init(&cal);
    cal.measurement_count = 15;
    EXPECT_EQ(15u, TRICAL_measurement_count_get(&cal));
}

/* The count is of valid measurements, whether or not they were needed */
TEST(TRICAL, MeasurementCountCountsValidMeasurements) {
    const float nan = std::numeric_limits<float>::quiet_NaN();
    TRICAL_instance_t cal;
    float measurement[3] = { 1.0f, 0.0f, 0.0f };
    float bad[3] = { 1.0f, nan, 0.0f };

    TRICAL_init(&cal);
    for (int i = 0; i < 200; i++) {
        TRICAL_estimate_update(&cal, measurement, nullptr);
    }
    EXPECT_EQ(200u, TRICAL_measurement_count_get(&cal));

    TRICAL_estimate_update(&cal, bad, nullptr);
    EXPECT_EQ(200u, TRICAL_measurement_count_get(&cal));
}

/*
Check that estimate update works with a series of measurements that should
result in the identity matrix with zero bias
*/
TEST(TRICAL, EstimateUpdateIdentity) {
    TRICAL_instance_t cal;

    TRICAL_init(&cal);

    float measurement[3];
    unsigned int i;

    for (i = 0; i < 100; i++) {
        measurement[0] = 1.0;
        measurement[1] = 0.0;
        measurement[2] = 0.0;
        TRICAL_estimate_update(&cal, measurement, measurement);

        measurement[0] = 0.0;
        measurement[1] = 1.0;
        measurement[2] = 0.0;
        TRICAL_estimate_update(&cal, measurement, measurement);

        measurement[0] = 0.0;
        measurement[1] = 0.0;
        measurement[2] = 1.0;
        TRICAL_estimate_update(&cal, measurement, measurement);

        measurement[0] = -1.0;
        measurement[1] = 0.0;
        measurement[2] = 0.0;
        TRICAL_estimate_update(&cal, measurement, measurement);

        measurement[0] = 0.0;
        measurement[1] = -1.0;
        measurement[2] = 0.0;
        TRICAL_estimate_update(&cal, measurement, measurement);

        measurement[0] = 0.0;
        measurement[1] = 0.0;
        measurement[2] = -1.0;
        TRICAL_estimate_update(&cal, measurement, measurement);
    }

    float bias_estimate[3], scale_estimate[9];
    TRICAL_estimate_get(&cal, bias_estimate, scale_estimate);
    for (i = 0; i < 3; i++) {
        EXPECT_NEAR(0.0, bias_estimate[i], 1e-5);
    }
    for (i = 0; i < 9; i++) {
        EXPECT_NEAR(0.0, scale_estimate[i], 1e-5);
    }
}

/*
Check that estimate update works with a series of measurements that should
result in the identity matrix with a bias in the positive X axis
*/
TEST(TRICAL, EstimateUpdateBias) {
    TRICAL_instance_t cal;

    TRICAL_init(&cal);

    float measurement[3], ref[3];
    unsigned int i;

    for (i = 0; i < 200; i++) {
        measurement[0] = 2.0;
        measurement[1] = 0.0;
        measurement[2] = 0.0;
        ref[0] = 1.0;
        ref[1] = 0.0;
        ref[2] = 0.0;
        TRICAL_estimate_update(&cal, measurement, ref);

        measurement[0] = 1.0;
        measurement[1] = 1.0;
        measurement[2] = 0.0;
        ref[0] = 0.0;
        ref[1] = 1.0;
        ref[2] = 0.0;
        TRICAL_estimate_update(&cal, measurement, ref);

        measurement[0] = 1.0;
        measurement[1] = 0.0;
        measurement[2] = 1.0;
        ref[0] = 0.0;
        ref[1] = 0.0;
        ref[2] = 1.0;
        TRICAL_estimate_update(&cal, measurement, ref);

        measurement[0] = 0.0;
        measurement[1] = 0.0;
        measurement[2] = 0.0;
        ref[0] = -1.0;
        ref[1] = 0.0;
        ref[2] = 0.0;
        TRICAL_estimate_update(&cal, measurement, ref);

        measurement[0] = 1.0;
        measurement[1] = -1.0;
        measurement[2] = 0.0;
        ref[0] = 0.0;
        ref[1] = -1.0;
        ref[2] = 0.0;
        TRICAL_estimate_update(&cal, measurement, ref);

        measurement[0] = 1.0;
        measurement[1] = 0.0;
        measurement[2] = -1.0;
        ref[0] = 0.0;
        ref[1] = 0.0;
        ref[2] = -1.0;
        TRICAL_estimate_update(&cal, measurement, ref);
    }

    float bias_estimate[3], scale_estimate[9];
    TRICAL_estimate_get(&cal, bias_estimate, scale_estimate);
    EXPECT_NEAR(1.0, bias_estimate[0], 1e-4);
    EXPECT_NEAR(0.0, bias_estimate[1], 1e-4);
    EXPECT_NEAR(0.0, bias_estimate[2], 1e-4);
    for (i = 0; i < 9; i++) {
        EXPECT_NEAR(0.0, scale_estimate[i], 1e-4);
    }
}

/* NULL and the measurement itself both mean "no reference" */
TEST(TRICAL, ReferenceIsOptional) {
    TRICAL_instance_t with_null, with_alias;
    Random rng(1);
    const Sensor sensor = Sensor::random(rng, 1.0, 0.3, 0.1, 0.01);
    Tumble motion(rng);

    TRICAL_init(&with_null);
    TRICAL_init(&with_alias);

    for (int i = 0; i < 500; i++) {
        const Vec3 reading = sensor.read(motion.step(0.03), rng);
        float measurement[3] = { static_cast<float>(reading[0]),
                                 static_cast<float>(reading[1]),
                                 static_cast<float>(reading[2]) };

        TRICAL_estimate_update(&with_null, measurement, nullptr);
        TRICAL_estimate_update(&with_alias, measurement, measurement);
    }

    for (int i = 0; i < TRICAL_STATE_DIM; i++) {
        EXPECT_EQ(with_null.state[i], with_alias.state[i]);
    }

    /* without a reference the scale estimate is symmetric */
    EXPECT_EQ(with_null.state[4], with_null.state[6]);
    EXPECT_EQ(with_null.state[5], with_null.state[9]);
    EXPECT_EQ(with_null.state[8], with_null.state[10]);
}

/* Check that state estimates can be accessed */
TEST(TRICAL, EstimateGet) {
    TRICAL_instance_t cal;

    TRICAL_init(&cal);
    cal.state[0] = 3.0;
    cal.state[1] = 2.0;
    cal.state[2] = 1.0;
    cal.state[3] = 1.0;
    cal.state[4] = 0.5;
    cal.state[5] = 0.1;
    cal.state[6] = 0.5;
    cal.state[7] = 1.0;
    cal.state[8] = 0.5;
    cal.state[9] = 0.1;
    cal.state[10] = 0.5;
    cal.state[11] = 1.0;

    float bias_estimate[3], scale_estimate[9];
    TRICAL_estimate_get(&cal, bias_estimate, scale_estimate);
    EXPECT_FLOAT_EQ(3.0, bias_estimate[0]);
    EXPECT_FLOAT_EQ(2.0, bias_estimate[1]);
    EXPECT_FLOAT_EQ(1.0, bias_estimate[2]);

    EXPECT_FLOAT_EQ(1.0, scale_estimate[0]);
    EXPECT_FLOAT_EQ(0.5, scale_estimate[1]);
    EXPECT_FLOAT_EQ(0.1, scale_estimate[2]);

    EXPECT_FLOAT_EQ(0.5, scale_estimate[3]);
    EXPECT_FLOAT_EQ(1.0, scale_estimate[4]);
    EXPECT_FLOAT_EQ(0.5, scale_estimate[5]);

    EXPECT_FLOAT_EQ(0.1, scale_estimate[6]);
    EXPECT_FLOAT_EQ(0.5, scale_estimate[7]);
    EXPECT_FLOAT_EQ(1.0, scale_estimate[8]);
}

/*
Check that state and variance estimates can be accessed. Before any
measurements the variances are those of the prior.
*/
TEST(TRICAL, EstimateGetExt) {
    TRICAL_instance_t cal;

    TRICAL_init(&cal);
    TRICAL_norm_set(&cal, 2.0f);

    cal.state[0] = 3.0;
    cal.state[1] = 2.0;
    cal.state[2] = 1.0;
    cal.state[3] = 1.0;
    cal.state[4] = 0.5;
    cal.state[5] = 0.1;
    cal.state[6] = 0.5;
    cal.state[7] = 1.0;
    cal.state[8] = 0.5;
    cal.state[9] = 0.1;
    cal.state[10] = 0.5;
    cal.state[11] = 1.0;

    float bias_estimate[3], scale_estimate[9], bias_estimate_variance[3],
          scale_estimate_variance[9];
    TRICAL_estimate_get_ext(&cal, bias_estimate, scale_estimate,
                            bias_estimate_variance, scale_estimate_variance);

    EXPECT_FLOAT_EQ(3.0, bias_estimate[0]);
    EXPECT_FLOAT_EQ(2.0, bias_estimate[1]);
    EXPECT_FLOAT_EQ(1.0, bias_estimate[2]);

    EXPECT_FLOAT_EQ(1.0, scale_estimate[0]);
    EXPECT_FLOAT_EQ(0.5, scale_estimate[1]);
    EXPECT_FLOAT_EQ(0.1, scale_estimate[2]);

    EXPECT_FLOAT_EQ(0.5, scale_estimate[3]);
    EXPECT_FLOAT_EQ(1.0, scale_estimate[4]);
    EXPECT_FLOAT_EQ(0.5, scale_estimate[5]);

    EXPECT_FLOAT_EQ(0.1, scale_estimate[6]);
    EXPECT_FLOAT_EQ(0.5, scale_estimate[7]);
    EXPECT_FLOAT_EQ(1.0, scale_estimate[8]);

    /*
    Before there are any measurements: a tenth of the size of the ellipsoid
    in the bias and a tenth in the scale factors, in the units that the
    scale estimate implies -- the mean square of the singular values of
    (I + D), which is a third of the sum of the squares of its elements.
    */
    const float units2 = (3 * 2.0f * 2.0f + 4 * 0.5f * 0.5f +
                          2 * 0.1f * 0.1f) / 3.0f;
    for (int i = 0; i < 3; i++) {
        /* (0.1 * field norm)^2, in the units of the measurements */
        EXPECT_NEAR(0.04f / units2, bias_estimate_variance[i], 1e-6f);
    }
    for (int i = 0; i < 9; i++) {
        EXPECT_NEAR(0.01f * units2, scale_estimate_variance[i], 1e-6f);
    }
}

/*
Check that measurement calibration can be performed based on the instance
state.
*/
TEST(TRICAL, Calibrate) {
    TRICAL_instance_t cal;

    TRICAL_init(&cal);
    cal.state[0] = 3.0;
    cal.state[1] = 2.0;
    cal.state[2] = 1.0;
    cal.state[3] = 1.0;
    cal.state[4] = 0.5;
    cal.state[5] = 0.1;
    cal.state[6] = 0.5;
    cal.state[7] = 1.0;
    cal.state[8] = 0.5;
    cal.state[9] = 0.1;
    cal.state[10] = 0.5;
    cal.state[11] = 1.0;

    float measurement[3] = { 1.0, 2.0, 3.0 }, result[3];
    TRICAL_measurement_calibrate(&cal, measurement, result);
    EXPECT_FLOAT_EQ(-3.8, result[0]);
    EXPECT_FLOAT_EQ(0.0, result[1]);
    EXPECT_FLOAT_EQ(3.8, result[2]);
}

/*
Measurement calibration: identity, bias only, scale only, and symmetric
scale, as (I + D)(B - b).
*/
TEST(TRICAL, CalibrateCases) {
    TRICAL_instance_t cal;
    float measurement[3] = { 1.0, 2.0, 3.0 }, result[3];

    TRICAL_init(&cal);
    TRICAL_measurement_calibrate(&cal, measurement, result);
    EXPECT_FLOAT_EQ(1.0, result[0]);
    EXPECT_FLOAT_EQ(2.0, result[1]);
    EXPECT_FLOAT_EQ(3.0, result[2]);

    /* bias only */
    cal.state[0] = 1.0;
    cal.state[1] = 2.0;
    cal.state[2] = 3.0;
    TRICAL_measurement_calibrate(&cal, measurement, result);
    EXPECT_FLOAT_EQ(0.0, result[0]);
    EXPECT_FLOAT_EQ(0.0, result[1]);
    EXPECT_FLOAT_EQ(0.0, result[2]);

    /* D = I, so the scale matrix is 2 I */
    TRICAL_init(&cal);
    cal.state[3] = cal.state[7] = cal.state[11] = 1.0;
    TRICAL_measurement_calibrate(&cal, measurement, result);
    EXPECT_FLOAT_EQ(2.0, result[0]);
    EXPECT_FLOAT_EQ(4.0, result[1]);
    EXPECT_FLOAT_EQ(6.0, result[2]);

    /* symmetric D */
    const float d[9] = { 1.0, 0.5, 0.1, 0.5, 1.0, 0.5, 0.1, 0.5, 1.0 };
    std::memcpy(&cal.state[3], d, sizeof(d));
    TRICAL_measurement_calibrate(&cal, measurement, result);
    EXPECT_FLOAT_EQ(3.3, result[0]);
    EXPECT_FLOAT_EQ(6.0, result[1]);
    EXPECT_FLOAT_EQ(7.1, result[2]);
}

/* The input and output of calibration may be the same vector */
TEST(TRICAL, CalibrateInPlace) {
    TRICAL_instance_t cal;

    TRICAL_init(&cal);
    const float state[12] = { 3.0, 2.0, 1.0, 1.0, 0.5, 0.1, 0.5, 1.0, 0.5,
                              0.1, 0.5, 1.0 };
    std::memcpy(cal.state, state, sizeof(state));

    float measurement[3] = { 1.0, 2.0, 3.0 };
    TRICAL_measurement_calibrate(&cal, measurement, measurement);
    EXPECT_FLOAT_EQ(-3.8, measurement[0]);
    EXPECT_FLOAT_EQ(0.0, measurement[1]);
    EXPECT_FLOAT_EQ(3.8, measurement[2]);
}

/* Reset discards what was learned, but not the configuration or count */
TEST(TRICAL, Reset) {
    TRICAL_instance_t cal;
    Random rng(2);
    const Sensor sensor = Sensor::random(rng, 3.0, 0.4, 0.1, 0.01);

    TRICAL_init(&cal);
    TRICAL_norm_set(&cal, 3.0f);
    TRICAL_noise_set(&cal, 0.05f);
    tumble(cal, sensor, rng, 2000);
    ASSERT_LT(support::errors_of(cal, sensor).magnitude, 5e-3);

    TRICAL_reset(&cal);
    EXPECT_FLOAT_EQ(3.0f, TRICAL_norm_get(&cal));
    EXPECT_FLOAT_EQ(0.05f, TRICAL_noise_get(&cal));
    EXPECT_EQ(2000u, TRICAL_measurement_count_get(&cal));
    for (int i = 0; i < TRICAL_STATE_DIM; i++) {
        EXPECT_EQ(0.0f, cal.state[i]);
    }

    /* and it learns again from scratch: one sample is not enough */
    support::update(cal, sensor.read(rng.direction(), rng));
    EXPECT_GT(support::errors_of(cal, sensor).magnitude, 0.1);
    tumble(cal, sensor, rng, 2000);
    EXPECT_LT(support::errors_of(cal, sensor).magnitude, 5e-3);
}

/* Init on a used instance restores the defaults */
TEST(TRICAL, InitResetsEverything) {
    TRICAL_instance_t cal, fresh;
    Random rng(3);
    const Sensor sensor = Sensor::random(rng, 1.0, 0.4, 0.1, 0.01);

    TRICAL_init(&cal);
    tumble(cal, sensor, rng, 500);
    TRICAL_init(&cal);

    TRICAL_init(&fresh);
    EXPECT_EQ(0, std::memcmp(&cal, &fresh, sizeof(cal)));
}

/* Setting the norm before any measurements does not touch the estimate */
TEST(TRICAL, NormSetBeforeMeasurements) {
    TRICAL_instance_t cal;

    TRICAL_init(&cal);
    TRICAL_norm_set(&cal, 60.0f);
    for (int i = 0; i < TRICAL_STATE_DIM; i++) {
        EXPECT_EQ(0.0f, cal.state[i]);
    }
}

/*
Changing the norm once the calibration is known rescales the output, and
leaves the bias alone
*/
TEST(TRICAL, NormSetRescalesEstimate) {
    TRICAL_instance_t cal;
    Random rng(4);
    const Sensor sensor = Sensor::random(rng, 1.0, 0.4, 0.1, 0.005);

    TRICAL_init(&cal);
    tumble(cal, sensor, rng, 3000);

    float bias_before[3], scale[9], bias_after[3];
    TRICAL_estimate_get(&cal, bias_before, scale);

    TRICAL_norm_set(&cal, 2.5f);
    TRICAL_estimate_get(&cal, bias_after, scale);
    for (int i = 0; i < 3; i++) {
        EXPECT_EQ(bias_before[i], bias_after[i]);
    }

    /* immediately... */
    Tumble motion(rng);
    for (int i = 0; i < 50; i++) {
        const Vec3 reading = sensor.read(motion.step(0.3), rng);
        float measurement[3] = { static_cast<float>(reading[0]),
                                 static_cast<float>(reading[1]),
                                 static_cast<float>(reading[2]) };
        float out[3];

        TRICAL_measurement_calibrate(&cal, measurement, out);
        EXPECT_NEAR(2.5, std::sqrt(out[0] * out[0] + out[1] * out[1] +
                                   out[2] * out[2]), 0.06);
    }

    /* ...and after further measurements */
    Sensor rescaled = sensor;
    rescaled.field_norm = 2.5;
    rescaled.noise = 0.005;
    for (int i = 0; i < 9; i++) {
        rescaled.a[i] *= 2.5;
    }
    tumble(cal, sensor, rng, 500);
    EXPECT_LT(support::errors_of(cal, rescaled).magnitude, 3e-3);
}

/*
A calibration written to `state` is used straight away, and is the starting
point for further estimation
*/
TEST(TRICAL, StateWriteRestartsFromWrittenCalibration) {
    TRICAL_instance_t cal;
    Random rng(5);
    Sensor sensor = Sensor::random(rng, 1.0, 0.8, 0.15, 0.01);
    sensor.misalign(rng.direction(), 0.2);

    TRICAL_init(&cal);

    /* something to be forgotten */
    {
        Random other(55);
        const Sensor different = Sensor::random(other, 1.0, 0.5, 0.1, 0.01);
        tumble(cal, different, other, 1000);
    }

    /* write the correct calibration, as a restore from storage would */
    for (int i = 0; i < 3; i++) {
        cal.state[i] = static_cast<float>(sensor.b[i]);
    }
    for (int i = 0; i < 9; i++) {
        cal.state[3 + i] = static_cast<float>(sensor.a[i]) -
                           ((i % 4 == 0) ? 1.0f : 0.0f);
    }
    EXPECT_LT(support::errors_of(cal, sensor).vector, 1e-6);

    /*
    A few measurements in one place cannot say much: the written calibration
    (including its rotation) must survive them.
    */
    const Vec3 direction = rng.direction();
    for (int i = 0; i < 300; i++) {
        support::update(cal, sensor.read(direction, rng));
    }
    EXPECT_LT(support::errors_of(cal, sensor).vector, 0.02);

    /* and it stays right as the data come in */
    tumble(cal, sensor, rng, 2000);
    EXPECT_LT(support::errors_of(cal, sensor).magnitude, 3e-3);
    EXPECT_LT(support::errors_of(cal, sensor).vector, 3e-3);
}

/*
The estimate is refreshed when the sensor has moved far enough to start a new
cell of measurements -- but not more often than every few measurements, once
there are plenty, however fast the sensor turns.
*/
TEST(TRICAL, RefitsAreSpacedOut) {
    TRICAL_instance_t cal;
    Random rng(9);
    const Sensor sensor = Sensor::random(rng, 1.0, 0.3, 0.1, 0.01);
    float previous[TRICAL_STATE_DIM];
    int first = 0, early = 0, late = 0;

    TRICAL_init(&cal);
    std::memcpy(previous, cal.state, sizeof(previous));

    /*
    Every reading is somewhere else entirely (and is used one update later,
    when the next reading has shown that it was not a glitch).
    */
    for (int i = 0; i < 1000; i++) {
        support::update(cal, sensor.read(rng.direction(), rng));
        if (std::memcmp(previous, cal.state, sizeof(previous)) != 0) {
            std::memcpy(previous, cal.state, sizeof(previous));
            (i < 9 ? first : (i < 33 ? early : late))++;
        }
    }

    /* there is nothing to fit until there are as many readings as unknowns */
    EXPECT_EQ(0, first);
    /* then every reading counts at once... */
    EXPECT_GE(early, 16);
    /* ...until there are plenty */
    EXPECT_LE(late, 967 / 8 + 24);
    EXPECT_GE(late, 967 / 8 - 24);
    EXPECT_LT(support::errors_of(cal, sensor).magnitude, 2e-3);
}

/*
The norm may be changed as often as you like: `state` follows it, and comes
back to exactly what it was when the norm does.
*/
TEST(TRICAL, NormChangesAreExactlyReversible) {
    TRICAL_instance_t cal;
    Random rng(10);
    const Sensor sensor = Sensor::random(rng, 1.0, 0.3, 0.1, 0.01);
    float first[TRICAL_STATE_DIM], doubled[TRICAL_STATE_DIM];

    TRICAL_init(&cal);
    tumble(cal, sensor, rng, 2000);

    TRICAL_norm_set(&cal, 2.0f);
    std::memcpy(doubled, cal.state, sizeof(doubled));
    TRICAL_norm_set(&cal, 1.0f);
    std::memcpy(first, cal.state, sizeof(first));
    EXPECT_LT(support::errors_of(cal, sensor).magnitude, 2e-3);

    /* the bias does not depend on the norm; the scale matrix follows it */
    for (int i = 0; i < 3; i++) {
        EXPECT_EQ(first[i], doubled[i]);
    }
    for (int i = 0; i < 9; i++) {
        const float identity = (i % 4 == 0) ? 1.0f : 0.0f;

        EXPECT_NEAR(2.0f * (first[3 + i] + identity),
                    doubled[3 + i] + identity, 1e-6f);
    }

    /* e.g. a field model that is evaluated afresh for every reading */
    for (int i = 0; i < 100000; i++) {
        TRICAL_norm_set(&cal, 1.0f + 0.3f * std::sin(0.001f * (float)i));
    }
    TRICAL_norm_set(&cal, 2.0f);
    EXPECT_EQ(0, std::memcmp(doubled, cal.state, sizeof(doubled)));
    TRICAL_norm_set(&cal, 1.0f);
    EXPECT_EQ(0, std::memcmp(first, cal.state, sizeof(first)));
}

/*
Likewise for a sensor in raw units, where the scale matrix is small and
`state` holds it with little precision to spare.
*/
TEST(TRICAL, NormChangesDoNotDegradeRawUnitCalibration) {
    TRICAL_instance_t cal;
    Random rng(11);
    Sensor sensor = Sensor::random(rng, 1.0, 0.0, 0.1, 0.0);

    for (int i = 0; i < 9; i++) {
        sensor.a[i] /= 1500.0;
    }
    sensor.b = Vec3{{ 420.0, -310.0, 480.0 }};
    sensor.noise = 10.0;

    TRICAL_init(&cal);
    cal.state[3] = cal.state[7] = cal.state[11] = 1.0f / 1500.0f - 1.0f;
    tumble(cal, sensor, rng, 3000);
    const double before = support::errors_of(cal, sensor).magnitude;

    for (int i = 0; i < 100000; i++) {
        TRICAL_norm_set(&cal, 1.0f + 0.0001f * (float)(i % 100));
    }
    TRICAL_norm_set(&cal, 1.0f);
    EXPECT_LT(support::errors_of(cal, sensor).magnitude, before + 2e-4);
}

/*
A sensor that reports in raw units: the initial scale estimate written to
`state` says roughly how they relate to the units of the field norm, and
everything follows from that.
*/
TEST(TRICAL, StateWriteOfUnitConversion) {
    TRICAL_instance_t cal;
    Random rng(8);

    /* about 1500 counts per unit of field, and a bias of 700 counts */
    Sensor sensor = Sensor::random(rng, 1.0, 0.0, 0.1, 0.0);
    for (int i = 0; i < 9; i++) {
        sensor.a[i] /= 1500.0;
    }
    sensor.b = Vec3{{ 420.0, -310.0, 480.0 }};
    sensor.noise = 10.0;

    TRICAL_init(&cal);
    cal.state[3] = cal.state[7] = cal.state[11] = 1.0f / 1500.0f - 1.0f;

    tumble(cal, sensor, rng, 3000);
    EXPECT_LT(support::errors_of(cal, sensor).magnitude, 2e-3);
    EXPECT_NEAR(10.0, TRICAL_noise_estimate_get(&cal), 1.5);

    /* at rest it costs nothing: few of these readings are used */
    const unsigned int used = static_cast<unsigned int>(cal.priv.moments[0]);
    const Vec3 direction = rng.direction();
    for (int i = 0; i < 5000; i++) {
        support::update(cal, sensor.read(direction, rng));
    }
    EXPECT_LE(static_cast<unsigned int>(cal.priv.moments[0]), used + 128);
}

/* Garbage written to `state` is discarded at the next update */
TEST(TRICAL, StateWriteOfGarbageIsIgnored) {
    TRICAL_instance_t cal;
    Random rng(6);
    const Sensor sensor = Sensor::random(rng, 1.0, 0.3, 0.1, 0.01);

    TRICAL_init(&cal);
    tumble(cal, sensor, rng, 2000);

    float good[TRICAL_STATE_DIM];
    std::memcpy(good, cal.state, sizeof(good));

    /*
    Not finite. A reading that does not prompt a refit (here, one that is not
    usable itself) leaves exactly the calibration that was there before.
    */
    const float unusable[3] = { 0.0f, 0.0f,
                                std::numeric_limits<float>::infinity() };
    cal.state[4] = std::numeric_limits<float>::quiet_NaN();
    TRICAL_estimate_update(&cal, unusable, nullptr);
    EXPECT_EQ(0, std::memcmp(good, cal.state, sizeof(good)));

    cal.state[0] = std::numeric_limits<float>::infinity();
    support::update(cal, sensor.read(rng.direction(), rng));
    EXPECT_TRUE(support::is_finite(support::calibration_of(cal)));
    EXPECT_LT(support::errors_of(cal, sensor).magnitude, 5e-3);

    /* singular scale matrix (I + D = 0) */
    cal.state[3] = cal.state[7] = cal.state[11] = -1.0f;
    cal.state[4] = cal.state[5] = cal.state[6] = 0.0f;
    cal.state[8] = cal.state[9] = cal.state[10] = 0.0f;
    support::update(cal, sensor.read(rng.direction(), rng));
    EXPECT_TRUE(support::is_finite(support::calibration_of(cal)));
    EXPECT_LT(support::errors_of(cal, sensor).magnitude, 5e-3);
}

/* Instances are plain data: a copy carries on exactly like the original */
TEST(TRICAL, InstanceCanBeCopied) {
    TRICAL_instance_t original, copy;
    Random rng(7);
    const Sensor sensor = Sensor::random(rng, 1.0, 0.3, 0.1, 0.01);

    TRICAL_init(&original);
    tumble(original, sensor, rng, 700);
    std::memcpy(&copy, &original, sizeof(copy));

    Tumble motion(rng);
    for (int i = 0; i < 700; i++) {
        const Vec3 direction = motion.step(0.03);
        const Vec3 reading = sensor.read(direction, rng);

        support::update(original, reading, &direction);
        support::update(copy, reading, &direction);
    }
    EXPECT_EQ(0, std::memcmp(&copy, &original, sizeof(copy)));
}
