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

/*
Robustness: whatever is thrown at it, the estimator must not produce
non-finite output, and must not let bad input damage a good calibration.
*/

#include <gtest/gtest.h>

#include <cfloat>
#include <climits>
#include <cmath>
#include <cstring>
#include <limits>

/* C++ has no 'restrict' qualifier. Just ignore it. */
#define restrict

#include "TRICAL.h"
#include "fit.h"
#include "support.h"

using support::Random;
using support::Sensor;
using support::Tumble;
using support::Vec3;

namespace {

const double kDegree = 0.017453292519943295;
const float kNaN = std::numeric_limits<float>::quiet_NaN();
const float kInf = std::numeric_limits<float>::infinity();

/* A calibrated instance to abuse, and the sensor it belongs to. */
struct Calibrated {
    TRICAL_instance_t cal;
    Sensor sensor;
    Random rng;
    Tumble motion;

    explicit Calibrated(uint32_t seed, double field_norm = 1.0) :
            sensor(), rng(seed), motion(rng) {
        sensor = Sensor::random(rng, field_norm, 0.4, 0.1, 0.01);
        TRICAL_init(&cal);
        TRICAL_norm_set(&cal, static_cast<float>(field_norm));
        tumble(3000);
    }

    void tumble(int count) {
        for (int i = 0; i < count; i++) {
            const Vec3 direction = motion.step(kDegree);
            support::update(cal, sensor.read(direction, rng), &direction);
        }
    }

    bool state_equals(const float state[TRICAL_STATE_DIM]) const {
        return std::memcmp(cal.state, state, sizeof(cal.state)) == 0;
    }
};

bool finite_state(const TRICAL_instance_t &cal) {
    for (int i = 0; i < TRICAL_STATE_DIM; i++) {
        if (!std::isfinite(cal.state[i])) {
            return false;
        }
    }
    return true;
}

/* The number of measurements that have been used. */
double used(const TRICAL_instance_t &cal) {
    return static_cast<double>(cal.priv.moments[0]);
}

/* The opposite direction: as far from `direction` as a reading can be. */
Vec3 opposite(const Vec3 &direction) {
    return Vec3{{ -direction[0], -direction[1], -direction[2] }};
}

} /* namespace */

/* Measurements with NaN or infinite components are ignored entirely */
TEST(Robustness, NonFiniteMeasurementsAreIgnored) {
    Calibrated c(41);
    TRICAL_instance_t before;
    const float bad[][3] = {
        { kNaN, 0.5f, 0.5f }, { 0.5f, kNaN, 0.5f }, { 0.5f, 0.5f, kNaN },
        { kInf, 0.5f, 0.5f }, { 0.5f, -kInf, 0.5f }, { kNaN, kNaN, kNaN },
        { kInf, kInf, kInf }
    };
    const float reference[3] = { 1.0f, 0.0f, 0.0f };

    std::memcpy(&before, &c.cal, sizeof(before));
    for (const float *measurement : bad) {
        TRICAL_estimate_update(&c.cal, measurement, nullptr);
        TRICAL_estimate_update(&c.cal, measurement, reference);
    }
    EXPECT_EQ(0, std::memcmp(&before, &c.cal, sizeof(before)));
}

/*
A reference that is not finite, or has no direction, is ignored -- but the
measurement that came with it is still used
*/
TEST(Robustness, UnusableReferencesAreIgnored) {
    Random rng(42);
    const Sensor sensor = Sensor::random(rng, 1.0, 0.4, 0.1, 0.01);
    TRICAL_instance_t with_bad, without;
    Tumble motion(rng);
    const float bad[][3] = {
        { kNaN, 0.0f, 0.0f }, { 0.0f, kInf, 0.0f }, { 0.0f, 0.0f, 0.0f },
        { 1e-30f, 0.0f, 0.0f }, { FLT_MAX, FLT_MAX, FLT_MAX }
    };

    TRICAL_init(&with_bad);
    TRICAL_init(&without);
    for (int i = 0; i < 2000; i++) {
        const Vec3 reading = sensor.read(motion.step(kDegree), rng);
        float measurement[3] = { static_cast<float>(reading[0]),
                                 static_cast<float>(reading[1]),
                                 static_cast<float>(reading[2]) };

        TRICAL_estimate_update(&with_bad, measurement, bad[i % 5]);
        TRICAL_estimate_update(&without, measurement, nullptr);
    }

    EXPECT_EQ(0, std::memcmp(with_bad.state, without.state,
                             sizeof(without.state)));
    EXPECT_LT(support::errors_of(with_bad, sensor).magnitude, 3e-3);
}

/* Absurd magnitudes neither overflow nor disturb the calibration */
TEST(Robustness, AbsurdMagnitudes) {
    Calibrated c(43);
    float good[TRICAL_STATE_DIM];
    const float values[] = { FLT_MAX, -FLT_MAX, 1e30f, -1e20f, 1e10f, 1e5f,
                             100.0f };

    std::memcpy(good, c.cal.state, sizeof(good));
    for (float value : values) {
        const float measurement[3] = { value, -value, value };
        const float reference[3] = { value, value, value };

        TRICAL_estimate_update(&c.cal, measurement, reference);
        EXPECT_TRUE(c.state_equals(good)) << value;
    }

    /* a fresh instance is no more easily upset */
    TRICAL_instance_t fresh;
    TRICAL_init(&fresh);
    for (float value : values) {
        const float measurement[3] = { value, -value, value };

        for (int i = 0; i < 400; i++) {
            TRICAL_estimate_update(&fresh, measurement, nullptr);
        }
        EXPECT_TRUE(finite_state(fresh)) << value;
    }
}

/* Tiny and zero measurements */
TEST(Robustness, TinyAndZeroMeasurements) {
    const float values[] = { 0.0f, -0.0f, FLT_MIN, 1e-40f, 1e-20f };

    for (float value : values) {
        TRICAL_instance_t cal;
        const float measurement[3] = { value, value, -value };

        TRICAL_init(&cal);
        for (int i = 0; i < 1000; i++) {
            TRICAL_estimate_update(&cal, measurement, measurement);
            ASSERT_TRUE(finite_state(cal)) << value;
        }

        float out[3];
        TRICAL_measurement_calibrate(&cal, measurement, out);
        EXPECT_TRUE(std::isfinite(out[0]) && std::isfinite(out[1]) &&
                    std::isfinite(out[2]));
    }
}

/* The same reading over and over */
TEST(Robustness, ConstantInput) {
    TRICAL_instance_t cal;
    const float measurement[3] = { 0.7f, -0.4f, 0.9f };

    TRICAL_init(&cal);
    for (int i = 0; i < 100000; i++) {
        TRICAL_estimate_update(&cal, measurement, nullptr);
    }
    ASSERT_TRUE(finite_state(cal));

    /* a single point says little; the estimate must stay near the start */
    for (int i = 0; i < TRICAL_STATE_DIM; i++) {
        EXPECT_LT(std::fabs(cal.state[i]), 0.3f) << "state " << i;
    }
}

/*
A calibrated sensor at rest: the estimate stays put, however long it rests.
(If the resting place is near the edge of the cell that the sensor was in
when it stopped, noise will sooner or later start one more cell there, and
with it refresh the estimate -- but only the once.)
*/
TEST(Robustness, RestDoesNotChangeEstimate) {
    const double noises[] = { 0.01, 0.05 };

    for (double noise : noises) {
        for (uint32_t seed = 440; seed < 446; seed++) {
            Random rng(seed);
            const Sensor sensor = Sensor::random(rng, 1.0, 0.3, 0.1, noise);
            Tumble motion(rng);
            TRICAL_instance_t cal;
            float good[TRICAL_STATE_DIM];

            TRICAL_init(&cal);
            for (int i = 0; i < 5000; i++) {
                support::update(cal, sensor.read(motion.step(kDegree), rng));
            }
            std::memcpy(good, cal.state, sizeof(good));
            const double used = static_cast<double>(cal.priv.moments[0]);

            const Vec3 direction = motion.direction();
            for (int i = 0; i < 200000; i++) {
                support::update(cal, sensor.read(direction, rng));
            }

            /* a small fraction of the statistical error of the estimate */
            for (int i = 0; i < TRICAL_STATE_DIM; i++) {
                EXPECT_NEAR(good[i], cal.state[i], 0.06 * noise)
                    << "noise " << noise << ", seed " << seed;
            }
            /* the rest of the current cell, and one more */
            EXPECT_LE(static_cast<double>(cal.priv.moments[0]),
                      used + 2.0 * cal.priv.cell_quota)
                << "noise " << noise << ", seed " << seed;
            EXPECT_LE(cal.priv.cell_quota, 64 * 5);
        }
    }
}

/* Once a cell has been started at the resting place, nothing changes */
TEST(Robustness, RestInOwnCellChangesNothing) {
    Calibrated c(44);
    float good[TRICAL_STATE_DIM];

    /* far enough from where the tumbling stopped to start a cell */
    Vec3 direction = c.motion.direction();
    for (int i = 0; i < 20; i++) {
        direction = c.motion.step(kDegree);
    }
    for (int i = 0; i < 10; i++) {
        support::update(c.cal, c.sensor.read(direction, c.rng));
    }
    std::memcpy(good, c.cal.state, sizeof(good));

    for (int i = 0; i < 200000; i++) {
        support::update(c.cal, c.sensor.read(direction, c.rng));
    }
    EXPECT_TRUE(c.state_equals(good));
}

/* A burst of wild readings is ignored */
TEST(Robustness, WildReadingsAreIgnored) {
    Calibrated c(45);
    float good[TRICAL_STATE_DIM];

    std::memcpy(good, c.cal.state, sizeof(good));
    for (int i = 0; i < 100; i++) {
        const Vec3 direction = c.rng.direction();
        const Vec3 wild = {{ 30.0 * direction[0], 30.0 * direction[1],
                             30.0 * direction[2] }};
        support::update(c.cal, wild);
    }
    EXPECT_TRUE(c.state_equals(good));

    /* and normal service resumes afterwards */
    c.tumble(500);
    EXPECT_LT(support::errors_of(c.cal, c.sensor).magnitude, 3e-3);
}

/* A few moderately bad readings among many good ones do limited damage */
TEST(Robustness, OccasionalDisturbancesDoLimitedDamage) {
    Calibrated c(46);

    for (int i = 0; i < 3000; i++) {
        const Vec3 direction = c.motion.step(kDegree);
        Vec3 reading = c.sensor.read(direction, c.rng);

        /* one reading in a hundred is off by 30% of the field */
        if (i % 100 == 0) {
            const Vec3 offset = c.rng.direction();
            for (int k = 0; k < 3; k++) {
                reading[k] += 0.3 * offset[k];
            }
        }
        support::update(c.cal, reading);
    }
    EXPECT_LT(support::errors_of(c.cal, c.sensor).magnitude, 1e-2);
}

/*
If the sensor's characteristics change completely, estimation starts again
by itself.
*/
TEST(Robustness, RecoversFromChangeOfSensor) {
    Calibrated c(47);
    Sensor changed = c.sensor;

    /* e.g. a strongly magnetised object now sits next to the sensor */
    changed.b = Vec3{{ c.sensor.b[0] + 15.0, c.sensor.b[1] - 20.0,
                       c.sensor.b[2] + 12.0 }};
    ASSERT_GT(support::errors_of(c.cal, changed).magnitude, 1.0);

    for (int i = 0; i < 4000; i++) {
        support::update(c.cal, changed.read(c.motion.step(kDegree), c.rng));
        ASSERT_TRUE(finite_state(c.cal));
    }
    EXPECT_LT(support::errors_of(c.cal, changed).magnitude, 3e-3);
}

/* Readings in units wildly different from the field norm are not used */
TEST(Robustness, WrongUnitsDoNoHarm) {
    TRICAL_instance_t cal;
    Random rng(48);
    Tumble motion(rng);

    /* norm left at 1.0, readings of magnitude 500 */
    TRICAL_init(&cal);
    for (int i = 0; i < 5000; i++) {
        const Vec3 d = motion.step(kDegree);
        const Vec3 reading = {{ 500.0 * d[0], 500.0 * d[1], 500.0 * d[2] }};

        support::update(cal, reading);
        ASSERT_TRUE(finite_state(cal));
    }
}

/*
An isolated reading that is not on the ellipsoid of the estimate is not
used, however near or far it is.
*/
TEST(Robustness, IsolatedOutliersAreIgnored) {
    const double distances[] = { 0.3, 0.6, 1.5, 2.0, 3.0, 5.0, 7.0 };
    Calibrated c(50);

    for (double distance : distances) {
        const Vec3 direction = c.rng.direction();
        const Vec3 outlier = {{ c.sensor.b[0] + distance * direction[0],
                                c.sensor.b[1] + distance * direction[1],
                                c.sensor.b[2] + distance * direction[2] }};
        const double before = used(c.cal);
        float good[TRICAL_STATE_DIM];

        std::memcpy(good, c.cal.state, sizeof(good));
        support::update(c.cal, outlier);
        EXPECT_TRUE(c.state_equals(good)) << "distance " << distance;
        EXPECT_EQ(before, used(c.cal)) << "distance " << distance;

        c.tumble(100);
    }
    EXPECT_LT(support::errors_of(c.cal, c.sensor).magnitude, 2e-3);
}

/*
A reading that is somewhere else when the next one is back where the sensor
was is a glitch, not a movement, and is not used -- even if it is a
perfectly plausible reading in itself. Two in a row are a movement.
*/
TEST(Robustness, GlitchIsNotUsedButMoveIs) {
    Calibrated c(51);
    const Vec3 here = c.motion.direction(), there = opposite(here);

    for (int i = 0; i < 500; i++) {
        support::update(c.cal, c.sensor.read(here, c.rng));
    }
    const double before = used(c.cal);

    /* there and back */
    support::update(c.cal, c.sensor.read(there, c.rng));
    EXPECT_EQ(before, used(c.cal));
    support::update(c.cal, c.sensor.read(here, c.rng));
    EXPECT_EQ(before, used(c.cal));

    /* there, and still there */
    support::update(c.cal, c.sensor.read(there, c.rng));
    EXPECT_EQ(before, used(c.cal));
    support::update(c.cal, c.sensor.read(there, c.rng));
    EXPECT_EQ(before + 2.0, used(c.cal));
}

/*
A sensor at rest that glitches now and then: the glitches are not used, and
neither do they make the readings at rest count all over again.
*/
TEST(Robustness, GlitchesAtRestDoNotAccumulate) {
    for (int kind = 0; kind < 4; kind++) {
        Calibrated c(52);
        const Vec3 here = c.motion.direction(), there = c.rng.direction();
        const double before = used(c.cal);

        for (int i = 1; i <= 100000; i++) {
            Vec3 reading = c.sensor.read(here, c.rng);

            if (i % 100 == 0) {
                switch (kind) {
                case 0:
                    /* nowhere near anything */
                    reading = Vec3{{ 30.0, -30.0, 30.0 }};
                    break;
                case 1:
                    /* not a number at all */
                    reading[1] = kNaN;
                    break;
                case 2:
                    /* a reading from somewhere else: the same place */
                    reading = c.sensor.read(there, c.rng);
                    break;
                default:
                    /* a reading from somewhere else: anywhere */
                    reading = c.sensor.read(c.rng.direction(), c.rng);
                    break;
                }
            }
            support::update(c.cal, reading);
        }

        /* at most what the cell at the resting place had left to take */
        EXPECT_LE(used(c.cal), before + c.cal.priv.cell_quota)
            << "kind " << kind;
        EXPECT_LT(support::errors_of(c.cal, c.sensor).magnitude, 2e-3)
            << "kind " << kind;
    }
}

/*
Going back and forth between a few places does not keep adding to the weight
of the readings from those places.
*/
TEST(Robustness, RevisitedPlacesDoNotAccumulate) {
    Calibrated c(53);
    const Vec3 places[3] = { c.motion.direction(), c.rng.direction(),
                             c.rng.direction() };
    const double before = used(c.cal);

    for (int cycle = 0; cycle < 300; cycle++) {
        for (const Vec3 &place : places) {
            for (int i = 0; i < 100; i++) {
                support::update(c.cal, c.sensor.read(place, c.rng));
            }
        }
    }
    EXPECT_LE(used(c.cal), before + 3.0 * c.cal.priv.cell_quota);
    EXPECT_LT(support::errors_of(c.cal, c.sensor).magnitude, 2e-3);
}

/*
A sensor that is stuck, saturated or disconnected for a while (the same
reading over and over) changes nothing, however long it lasts.
*/
TEST(Robustness, StuckValuesChangeNothing) {
    const float values[] = { 32.767f, 0.0f, -3.0e5f, 0.25f };
    Calibrated c(54);
    const double before = used(c.cal);
    float good[TRICAL_STATE_DIM];

    std::memcpy(good, c.cal.state, sizeof(good));
    for (float value : values) {
        const float stuck[3] = { value, value, value };

        for (int i = 0; i < 2000; i++) {
            TRICAL_estimate_update(&c.cal, stuck, nullptr);
        }
        EXPECT_TRUE(c.state_equals(good)) << "value " << value;
        EXPECT_EQ(before, used(c.cal)) << "value " << value;
    }

    c.tumble(500);
    EXPECT_LT(support::errors_of(c.cal, c.sensor).magnitude, 2e-3);
}

/*
Nor does a long burst of readings from all over the place: they are not
from a sensor whose calibration could be learned, so what has been learned
is kept.
*/
TEST(Robustness, GarbageBurstChangesNothing) {
    Calibrated c(55);
    const double before = used(c.cal);
    float good[TRICAL_STATE_DIM];

    std::memcpy(good, c.cal.state, sizeof(good));
    for (int i = 0; i < 5000; i++) {
        const Vec3 direction = c.rng.direction();
        const double magnitude = 20.0 + 40.0 * c.rng.uniform();

        support::update(c.cal, Vec3{{ magnitude * direction[0],
                                      magnitude * direction[1],
                                      magnitude * direction[2] }});
    }
    EXPECT_TRUE(c.state_equals(good));
    EXPECT_EQ(before, used(c.cal));

    c.tumble(500);
    EXPECT_GT(used(c.cal), before);
    EXPECT_LT(support::errors_of(c.cal, c.sensor).magnitude, 2e-3);
}

/*
A calibration written to `state` is left exactly as it is until the sensor
has moved enough for there to be something to learn.
*/
TEST(Robustness, RestoredCalibrationIsNotDisturbedAtRest) {
    Random rng(56);
    const Sensor sensor = Sensor::random(rng, 1.0, 0.4, 0.1, 0.01);
    const Vec3 direction = rng.direction();
    TRICAL_instance_t cal;
    float written[TRICAL_STATE_DIM];

    TRICAL_init(&cal);
    for (int i = 0; i < 3; i++) {
        cal.state[i] = static_cast<float>(sensor.b[i]);
    }
    for (int i = 0; i < 9; i++) {
        cal.state[3 + i] = static_cast<float>(sensor.a[i]) -
                           ((i % 4 == 0) ? 1.0f : 0.0f);
    }
    std::memcpy(written, cal.state, sizeof(written));

    for (int i = 0; i < 1000; i++) {
        support::update(cal, sensor.read(direction, rng));
        ASSERT_EQ(0, std::memcmp(written, cal.state, sizeof(written)))
            << "reading " << i;
    }
}

/* A field norm that makes no sense, however it got there, does no harm */
TEST(Robustness, UnusableNormIsInert) {
    const float norms[] = { 0.0f, -1.0f, kNaN, kInf, 1e30f, 1e-30f };

    for (float norm : norms) {
        Random rng(57);
        const Sensor sensor = Sensor::random(rng, 1.0, 0.3, 0.1, 0.01);
        Tumble motion(rng);
        TRICAL_instance_t cal;

        TRICAL_init(&cal);
        cal.field_norm = norm;
        for (int i = 0; i < 500; i++) {
            float bias[3], scale[9], bias_variance[3], scale_variance[9];

            support::update(cal, sensor.read(motion.step(kDegree), rng));
            ASSERT_TRUE(finite_state(cal)) << "norm " << norm;

            TRICAL_estimate_get_ext(&cal, bias, scale, bias_variance,
                                    scale_variance);
            for (int k = 0; k < 3; k++) {
                ASSERT_TRUE(std::isfinite(bias_variance[k]))
                    << "norm " << norm;
            }
            for (int k = 0; k < 9; k++) {
                ASSERT_TRUE(std::isfinite(scale_variance[k]))
                    << "norm " << norm;
            }
        }
        for (int i = 0; i < TRICAL_STATE_DIM; i++) {
            EXPECT_EQ(0.0f, cal.state[i]) << "norm " << norm;
        }
        EXPECT_EQ(0.0, used(cal)) << "norm " << norm;

        /* and a sensible one puts things right */
        cal.field_norm = 1.0f;
        for (int i = 0; i < 3000; i++) {
            support::update(cal, sensor.read(motion.step(kDegree), rng));
        }
        EXPECT_LT(support::errors_of(cal, sensor).magnitude, 3e-3)
            << "norm " << norm;

        /* the same when estimation is under way */
        float good[TRICAL_STATE_DIM];
        const double before = used(cal);

        std::memcpy(good, cal.state, sizeof(good));
        cal.field_norm = norm;
        for (int i = 0; i < 500; i++) {
            float bias[3], scale[9], bias_variance[3], scale_variance[9];

            support::update(cal, sensor.read(motion.step(kDegree), rng));
            TRICAL_estimate_get_ext(&cal, bias, scale, bias_variance,
                                    scale_variance);
            for (int k = 0; k < 3; k++) {
                ASSERT_TRUE(std::isfinite(bias_variance[k]))
                    << "norm " << norm;
            }
            for (int k = 0; k < 9; k++) {
                ASSERT_TRUE(std::isfinite(scale_variance[k]))
                    << "norm " << norm;
            }
        }
        EXPECT_EQ(0, std::memcmp(good, cal.state, sizeof(good)))
            << "norm " << norm;
        EXPECT_EQ(before, used(cal)) << "norm " << norm;

        cal.field_norm = 1.0f;
        for (int i = 0; i < 500; i++) {
            support::update(cal, sensor.read(motion.step(kDegree), rng));
        }
        EXPECT_GT(used(cal), before) << "norm " << norm;
        EXPECT_LT(support::errors_of(cal, sensor).magnitude, 3e-3)
            << "norm " << norm;
    }
}

/*
An initial scale estimate for measurements so large or so small that the
squares of distances between them could not be represented is refused, and
such measurements do no harm: they are too far away to be used, or too
close together to be told apart.
*/
TEST(Robustness, UnrepresentableUnitsAreRefused) {
    const float sizes[] = { 1e-25f, 1e-17f, 1e17f, 1e25f, 1e36f };

    for (float size : sizes) {
        Random rng(60);
        Tumble motion(rng);
        TRICAL_instance_t cal;

        /* readings of about this size, and the conversion to match */
        TRICAL_init(&cal);
        cal.state[3] = cal.state[7] = cal.state[11] = 1.0f / size - 1.0f;
        for (int i = 0; i < 2000; i++) {
            const Vec3 d = motion.step(kDegree);

            support::update(cal, Vec3{{ size * d[0], size * d[1],
                                        size * d[2] }});
            for (int k = 0; k < TRICAL_STATE_DIM; k++) {
                ASSERT_EQ(0.0f, cal.state[k]) << "size " << size;
            }
        }
        EXPECT_LE(used(cal), 1.0 * cal.priv.cell_quota) << "size " << size;
    }
}

/* The counters stop at their largest value; they do not start again */
TEST(Robustness, CountersSaturate) {
    Calibrated c(58);

    c.cal.measurement_count = UINT_MAX - 5u;
    c.cal.priv.reference_count = UINT_MAX - 5u;
    c.tumble(200);
    EXPECT_EQ(UINT_MAX, c.cal.measurement_count);
    EXPECT_EQ(UINT_MAX, c.cal.priv.reference_count);
    EXPECT_TRUE(finite_state(c.cal));
}

/*
A very long run. With single-precision sums the weight of the past is halved
now and then, which keeps the sums within the precision available; with
double precision there is no need.
*/
TEST(Robustness, VeryLongRun) {
    const bool wide_is_float = sizeof(trical_wide_t) == sizeof(float);
    Random rng(59);
    const Sensor sensor = Sensor::random(rng, 1.0, 0.3, 0.1, 0.01);
    Tumble motion(rng);
    TRICAL_instance_t cal;

    TRICAL_init(&cal);
    for (int i = 0; i < 400000; i++) {
        /* fast enough that every reading is used */
        support::update(cal, sensor.read(motion.step(10.0 * kDegree), rng));
    }
    ASSERT_TRUE(finite_state(cal));

    if (wide_is_float) {
        EXPECT_LE(used(cal), static_cast<double>(TRICAL_COUNT_MAX));
        EXPECT_GE(used(cal), static_cast<double>(TRICAL_COUNT_MAX) / 2.0);
    } else {
        EXPECT_GT(used(cal), 390000.0);
    }
    EXPECT_LT(support::errors_of(cal, sensor).magnitude,
              wide_is_float ? 5e-4 : 1e-4);
    /* (in single precision the residual of the fit is itself uncertain) */
    EXPECT_NEAR(0.01, TRICAL_noise_estimate_get(&cal),
                wide_is_float ? 1.5e-3 : 5e-4);
}

/* Random bit patterns as input: nothing non-finite ever comes out */
TEST(Robustness, Fuzz) {
    std::mt19937 gen(49);
    TRICAL_instance_t cal;

    TRICAL_init(&cal);
    for (int i = 0; i < 300000; i++) {
        float measurement[3], reference[3];
        uint32_t bits[6];

        for (int k = 0; k < 6; k++) {
            bits[k] = gen();
            /* a fair share of ordinary values among the garbage */
            if (gen() % 4 == 0) {
                const float ordinary = static_cast<float>(
                    static_cast<int32_t>(gen() % 4001) - 2000) / 1000.0f;
                std::memcpy(&bits[k], &ordinary, sizeof(ordinary));
            }
        }
        std::memcpy(measurement, &bits[0], sizeof(measurement));
        std::memcpy(reference, &bits[3], sizeof(reference));

        switch (gen() % 64) {
        case 0:
            TRICAL_reset(&cal);
            break;
        case 1:
            /* the caller scribbles on the state */
            std::memcpy(cal.state, &bits[0], 3 * sizeof(float));
            break;
        case 2:
            TRICAL_norm_set(&cal, 0.5f + static_cast<float>(gen() % 100));
            break;
        default:
            break;
        }

        TRICAL_estimate_update(&cal, measurement,
                               (gen() % 3 == 0) ? nullptr : reference);
        ASSERT_TRUE(finite_state(cal)) << "iteration " << i;

        if (i % 1000 == 0) {
            float bias[3], scale[9], bias_variance[3], scale_variance[9];

            TRICAL_estimate_get_ext(&cal, bias, scale, bias_variance,
                                    scale_variance);
            for (int k = 0; k < 3; k++) {
                ASSERT_TRUE(std::isfinite(bias_variance[k]));
                ASSERT_GE(bias_variance[k], 0.0f);
            }
            for (int k = 0; k < 9; k++) {
                ASSERT_TRUE(std::isfinite(scale_variance[k]));
                ASSERT_GE(scale_variance[k], 0.0f);
            }
        }
    }
}

/* After any amount of abuse, a reset gives a working estimator */
TEST(Robustness, ResetAfterGarbage) {
    Random rng(50);
    const Sensor sensor = Sensor::random(rng, 1.0, 0.4, 0.1, 0.01);
    TRICAL_instance_t cal;
    Tumble motion(rng);

    TRICAL_init(&cal);
    for (int i = 0; i < 5000; i++) {
        const float garbage[3] = { static_cast<float>(3.0 * rng.normal()),
                                   static_cast<float>(3.0 * rng.normal()),
                                   static_cast<float>(3.0 * rng.normal()) };
        TRICAL_estimate_update(&cal, garbage, nullptr);
    }
    ASSERT_TRUE(finite_state(cal));

    TRICAL_reset(&cal);
    for (int i = 0; i < 3000; i++) {
        support::update(cal, sensor.read(motion.step(kDegree), rng));
    }
    EXPECT_LT(support::errors_of(cal, sensor).magnitude, 3e-3);
}
