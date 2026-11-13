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
Statistical tests: how good are the estimates, measured on a simulated sensor
with known characteristics? The random numbers are reproducible, so these
tests are deterministic; the thresholds leave a margin of a factor of two or
more over the values observed.
*/

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

#include "TRICAL.h"
#include "support.h"

using support::Errors;
using support::Mat3;
using support::Random;
using support::Sensor;
using support::Tumble;
using support::Vec3;

namespace {

const bool kWideIsFloat = sizeof(TRICAL_WIDE_TYPE) == sizeof(float);

/* One degree, the rotation per sample of the "tumbling" sensor below. */
const double kDegree = 0.017453292519943295;

struct Scenario {
    TRICAL_instance_t cal;
    Sensor sensor;
    Random rng;
    Tumble motion;
    bool with_reference;
    double reference_error;
    std::vector<Vec3> readings;
    std::vector<Vec3> directions;

    Scenario(uint32_t seed, double field_norm, double bias,
            double scale_error, double noise) :
            sensor(), rng(seed), motion(rng), with_reference(false),
            reference_error(0.0) {
        sensor = Sensor::random(rng, field_norm, bias, scale_error, noise);
        TRICAL_init(&cal);
        TRICAL_norm_set(&cal, static_cast<float>(field_norm));
    }

    /* Starts from the correct calibration, as if restored from storage. */
    void restore() {
        for (int i = 0; i < 3; i++) {
            cal.state[i] = static_cast<float>(sensor.b[i]);
        }
        for (int i = 0; i < 9; i++) {
            cal.state[3 + i] = static_cast<float>(sensor.a[i]) -
                               ((i % 4 == 0) ? 1.0f : 0.0f);
        }
    }

    void feed(const Vec3 &direction) {
        const Vec3 reading = sensor.read(direction, rng);

        readings.push_back(reading);
        directions.push_back(direction);
        if (with_reference) {
            Vec3 reference = direction;
            if (reference_error > 0.0) {
                reference = support::mul(support::rotation(rng.direction(),
                    reference_error * rng.normal()), direction);
            }
            support::update(cal, reading, &reference);
        } else {
            support::update(cal, reading);
        }
    }

    /* tumbling about a wandering axis */
    void tumble(int count, double step = kDegree) {
        for (int i = 0; i < count; i++) {
            feed(motion.step(step));
        }
    }

    /* turning about one fixed axis */
    void turn(int count, const Vec3 &axis, double step = kDegree) {
        for (int i = 0; i < count; i++) {
            feed(motion.step_about(axis, step));
        }
    }

    /* not moving */
    void rest(int count) {
        const Vec3 direction = motion.direction();
        for (int i = 0; i < count; i++) {
            feed(direction);
        }
    }

    /* independent random orientations */
    void scatter(int count) {
        for (int i = 0; i < count; i++) {
            feed(rng.direction());
        }
    }

    Errors errors() const {
        return support::errors_of(cal, sensor);
    }

    /*
    RMS error in the magnitude of the calibrated field, relative to the
    field norm, at the orientations the sensor has been in (rather than over
    all orientations, as for errors()).
    */
    double visited_error() {
        Sensor noiseless = sensor;
        double sum = 0.0;

        noiseless.noise = 0.0;
        for (const Vec3 &direction : directions) {
            const Vec3 reading = noiseless.read(direction, rng);
            float measurement[3] = { static_cast<float>(reading[0]),
                                     static_cast<float>(reading[1]),
                                     static_cast<float>(reading[2]) };
            float out[3];

            TRICAL_measurement_calibrate(&cal, measurement, out);

            const double norm = std::sqrt(
                static_cast<double>(out[0]) * out[0] +
                static_cast<double>(out[1]) * out[1] +
                static_cast<double>(out[2]) * out[2]) / sensor.field_norm;
            sum += (norm - 1.0) * (norm - 1.0);
        }

        return std::sqrt(sum / static_cast<double>(directions.size()));
    }
};

} /* namespace */

/*
With measurements at independent random orientations, the estimate is as
accurate as a batch maximum-likelihood fit to the same data.
*/
TEST(Accuracy, MatchesMaximumLikelihood) {
    struct Case {
        double noise;
        int count;
        double limit;
    };
    const Case cases[] = {
        { 0.01, 100, 1.15 }, { 0.01, 1000, 1.15 }, { 0.001, 1000, 1.15 },
        { 0.05, 1000, 1.15 }
    };
    const int trials = 40;

    for (const Case &c : cases) {
        double sum_trical = 0.0, sum_ml = 0.0;

        for (int trial = 0; trial < trials; trial++) {
            Scenario s(1000 + trial, 1.0, 0.3, 0.1, c.noise);

            s.scatter(c.count);

            const Errors trical = s.errors();
            const Errors ml = support::errors_of(
                support::maximum_likelihood(s.readings, 1.0,
                                            support::truth_of(s.sensor)),
                s.sensor);

            sum_trical += trical.magnitude * trical.magnitude;
            sum_ml += ml.magnitude * ml.magnitude;
        }

        const double ratio = std::sqrt(sum_trical / sum_ml);
        RecordProperty("ratio", std::to_string(ratio));
        EXPECT_LT(ratio, kWideIsFloat ? 2.0 : c.limit)
            << "noise " << c.noise << ", " << c.count << " measurements";
        /* sanity: the yardstick itself is where theory says it should be */
        EXPECT_LT(std::sqrt(sum_ml / trials),
                  4.0 * c.noise / std::sqrt(c.count / 9.0));
    }
}

/* The usual case: a sensor tumbling smoothly, a degree per sample */
TEST(Accuracy, SmoothTumble) {
    for (int trial = 0; trial < 10; trial++) {
        Scenario s(2000 + trial, 1.0, 0.3, 0.1, 0.01);

        s.tumble(3000);
        EXPECT_LT(s.errors().magnitude, 2.5e-3) << "trial " << trial;
    }
}

/* Exact data, and the default settings: converges to working precision */
TEST(Accuracy, NoiseFree) {
    for (int trial = 0; trial < 10; trial++) {
        Scenario s(2100 + trial, 1.0, 0.3, 0.1, 0.0);

        s.tumble(3000);
        EXPECT_LT(s.errors().magnitude, kWideIsFloat ? 2e-3 : 2e-5)
            << "trial " << trial;
    }
}

/* Measurements in physical units (here microtesla) work just the same */
TEST(Accuracy, PhysicalUnits) {
    for (int trial = 0; trial < 10; trial++) {
        Scenario unit(2200 + trial, 1.0, 0.5, 0.1, 0.01);
        Scenario microtesla(2200 + trial, 50.0, 0.5, 0.1, 0.01);

        unit.tumble(3000);
        microtesla.tumble(3000);
        EXPECT_LT(microtesla.errors().magnitude, 2.5e-3) << "trial " << trial;
        EXPECT_NEAR(unit.errors().magnitude, microtesla.errors().magnitude,
                    1e-4);
    }
}

/* The noise setting only matters at the very start */
TEST(Accuracy, NoiseSettingIsNotCritical) {
    const float settings[] = { 1e-6f, 1e-3f, 1e-2f, 0.2f, 5.0f };

    for (int trial = 0; trial < 5; trial++) {
        for (float setting : settings) {
            Scenario s(2300 + trial, 1.0, 0.3, 0.1, 0.01);

            TRICAL_noise_set(&s.cal, setting);
            s.tumble(3000);
            EXPECT_LT(s.errors().magnitude, 2.5e-3)
                << "trial " << trial << ", noise setting " << setting;
        }
    }
}

/* The measurement noise is estimated from the data */
TEST(Accuracy, NoiseIsEstimated) {
    const double levels[] = { 0.001, 0.01, 0.03 };

    for (double level : levels) {
        /*
        With single-precision sums, rounding looks like noise of about this
        much, so less than that cannot be resolved.
        */
        if (kWideIsFloat && level < 0.005) {
            continue;
        }

        for (int trial = 0; trial < 5; trial++) {
            Scenario s(2400 + trial, 2.0, 0.3, 0.1, level);

            s.tumble(3000);
            EXPECT_NEAR(level * 2.0, TRICAL_noise_estimate_get(&s.cal),
                        0.15 * level * 2.0)
                << "noise " << level << ", trial " << trial;
        }
    }
}

/* A bias several times the size of the field, and strong scale errors */
TEST(Accuracy, LargeBias) {
    for (int trial = 0; trial < 10; trial++) {
        Scenario s(2500 + trial, 1.0, 5.0, 0.25, 0.01);

        s.tumble(5000);
        EXPECT_LT(s.errors().magnitude, 3e-3) << "trial " << trial;
    }
}

/* A sensor that starts at rest learns nothing wrong while it waits */
TEST(Accuracy, RestThenTumble) {
    for (int trial = 0; trial < 10; trial++) {
        Scenario s(2600 + trial, 1.0, 0.3, 0.1, 0.01);
        const double initial = s.errors().magnitude;

        s.rest(10000);
        EXPECT_LT(s.errors().magnitude, 1.5 * initial) << "trial " << trial;

        s.tumble(3000);
        EXPECT_LT(s.errors().magnitude, 2.5e-3) << "trial " << trial;
    }
}

/* Nor does a calibrated sensor lose anything while at rest */
TEST(Accuracy, TumbleThenRest) {
    for (int trial = 0; trial < 10; trial++) {
        Scenario s(2700 + trial, 1.0, 0.3, 0.1, 0.01);

        s.tumble(3000);
        const double before = s.errors().magnitude;

        s.rest(100000);
        EXPECT_LT(s.errors().magnitude, before + 1e-4) << "trial " << trial;
    }
}

/*
Turning about a single axis leaves part of the calibration undetermined: the
estimate must stay sensible while that lasts, and come right once the sensor
moves freely.
*/
TEST(Accuracy, TurnThenTumble) {
    const Vec3 axis = {{ 0.2, 0.1, 1.0 }};

    for (int trial = 0; trial < 10; trial++) {
        Scenario s(2800 + trial, 1.0, 0.3, 0.1, 0.01);
        const double initial = s.errors().magnitude;

        s.turn(3000, axis);
        EXPECT_LT(s.errors().magnitude, 2.0 * initial) << "trial " << trial;

        s.tumble(3000);
        EXPECT_LT(s.errors().magnitude, 2.5e-3) << "trial " << trial;
    }
}

/* A long spell of turning about one axis does not undo a good calibration */
TEST(Accuracy, TumbleThenLongTurn) {
    const Vec3 axis = {{ 0.2, 0.1, 1.0 }};

    for (int trial = 0; trial < 4; trial++) {
        Scenario s(2900 + trial, 1.0, 0.3, 0.1, 0.02);

        s.tumble(1500);
        const double before = s.errors().magnitude;

        s.turn(60000, axis);
        EXPECT_LT(s.errors().magnitude, std::max(2.0 * before, 6e-3))
            << "trial " << trial;
    }
}

/*
Turning about one axis only cannot determine the whole calibration. What it
does determine is found: the calibrated field has the right magnitude at the
orientations the sensor has been in. What it does not is left near the
initial estimate.
*/
TEST(Accuracy, TurnOnly) {
    const Vec3 axis = {{ 0.2, 0.1, 1.0 }};

    for (int trial = 0; trial < 10; trial++) {
        Scenario s(3600 + trial, 1.0, 0.3, 0.1, 0.01);

        s.turn(4000, axis);
        EXPECT_LT(s.visited_error(), 2e-3) << "trial " << trial;

        /* the initial estimate is in error by about 0.3 */
        EXPECT_LT(s.errors().magnitude, 0.35) << "trial " << trial;
    }
}

/*
If the initial estimate is right, it stays right: the measurements do not
contradict it, so what they do not determine is held firmly. (Noise would
otherwise push that part of the estimate about by a tenth of the field norm
and more, the measurements being unable to say that it is wrong.)
*/
TEST(Accuracy, TurnOnlyKeepsGoodInitialEstimate) {
    const Vec3 axis = {{ 0.2, 0.1, 1.0 }};
    const double noises[] = { 0.002, 0.01, 0.03 };

    for (double noise : noises) {
        for (int trial = 0; trial < 6; trial++) {
            Scenario s(4100 + trial, 1.0, 0.3, 0.1, noise);

            s.restore();
            for (int stage = 0; stage < 4; stage++) {
                s.turn(5000, axis);
                EXPECT_LT(s.errors().magnitude, 0.025)
                    << "noise " << noise << ", trial " << trial;
            }

            /* and it is no obstacle once the sensor moves freely */
            s.tumble(3000);
            EXPECT_LT(s.errors().magnitude, 0.3 * noise)
                << "noise " << noise << ", trial " << trial;
        }
    }
}

/* The same for readings from a small part of the ellipsoid */
TEST(Accuracy, SmallCoverageKeepsGoodInitialEstimate) {
    const double noises[] = { 0.0005, 0.002, 0.01, 0.03 };

    for (double noise : noises) {
        for (int trial = 0; trial < 6; trial++) {
            Scenario s(4200 + trial, 1.0, 0.3, 0.1, noise);
            const Vec3 centre = s.rng.direction();
            Vec3 offset = {{ 0.0, 0.0, 0.0 }};

            s.restore();

            /* wandering within about 20 degrees of one direction */
            for (int i = 0; i < 5000; i++) {
                Vec3 direction;
                double norm = 0.0;

                for (int k = 0; k < 3; k++) {
                    offset[k] += 0.02 * s.rng.normal();
                    offset[k] = std::max(-0.25, std::min(0.25, offset[k]));
                    direction[k] = centre[k] + offset[k];
                    norm += direction[k] * direction[k];
                }
                for (int k = 0; k < 3; k++) {
                    direction[k] /= std::sqrt(norm);
                }
                s.feed(direction);
            }
            EXPECT_LT(s.errors().magnitude, 0.03)
                << "noise " << noise << ", trial " << trial;
        }
    }
}

/* The same, with a bias several times the size of the field */
TEST(Accuracy, TurnOnlyWithLargeBias) {
    const Vec3 axis = {{ 0.2, 0.1, 1.0 }};

    for (int trial = 0; trial < 10; trial++) {
        Scenario s(3700 + trial, 50.0, 2.5, 0.1, 0.02);

        s.turn(4000, axis);
        EXPECT_LT(s.visited_error(), 6e-3) << "trial " << trial;

        /*
        The scale factors stay near their initial values, so the estimate
        is an ellipsoid of about the right size through the measurements:
        its centre cannot be further than twice that size from the truth.
        */
        EXPECT_LT(s.errors().magnitude, 1.5) << "trial " << trial;
    }
}

/* Turns about two different axes leave less undetermined */
TEST(Accuracy, TurnsAboutTwoAxes) {
    const Vec3 first = {{ 0.2, 0.1, 1.0 }}, second = {{ 1.0, 0.3, -0.1 }};

    for (int trial = 0; trial < 10; trial++) {
        Scenario s(3900 + trial, 1.0, 0.3, 0.1, 0.01);

        s.turn(1500, first);
        s.turn(1500, second);
        EXPECT_LT(s.visited_error(), 2e-3) << "trial " << trial;
        EXPECT_LT(s.errors().magnitude, 0.2) << "trial " << trial;
    }
}

/*
A handful of noisy measurements over a short arc say next to nothing about
the calibration, and must not be taken to: the estimate stays in the
neighbourhood of the initial one. (An unregularised fit to such data gives
a small ellipsoid hugging the arc, or worse.)
*/
TEST(Accuracy, ShortNoisyArc) {
    const Vec3 axis = {{ 0.2, 0.1, 1.0 }};

    for (int trial = 0; trial < 50; trial++) {
        Scenario s(3800 + trial, 1.0, 0.1, 0.03, 0.05);

        s.turn(20, axis);
        for (int i = 0; i < 3; i++) {
            EXPECT_LT(std::fabs(s.cal.state[i] - s.sensor.b[i]), 0.8)
                << "trial " << trial;
        }
        for (int i = 3; i < TRICAL_STATE_DIM; i++) {
            EXPECT_LT(std::fabs(s.cal.state[i]), 0.4) << "trial " << trial;
        }
    }
}

/*
A sensor whose readings are of a quite different size from the field norm is
learned without being told, so long as the readings are accepted at all (see
TRICAL_RANGE) and are far enough apart to be told from one another.
*/
TEST(Accuracy, ScaleFarFromInitialEstimate) {
    const double sizes[] = { 0.25, 0.4, 2.5, 5.0, 7.0 };

    for (double size : sizes) {
        for (int trial = 0; trial < 4; trial++) {
            Scenario s(4300 + trial, 1.0, 0.0, 0.1, 0.0);

            /* readings `size` times the field norm, noise 1% of that */
            for (int i = 0; i < 9; i++) {
                s.sensor.a[i] /= size;
            }
            s.sensor.b = Vec3{{ 0.3 * size, -0.2 * size, 0.1 * size }};
            s.sensor.noise = 0.01 * size;

            s.tumble(5000);
            EXPECT_LT(s.errors().magnitude, 2e-3)
                << "size " << size << ", trial " << trial;
        }
    }
}

/* A slowly turning sensor: many samples per degree */
TEST(Accuracy, SlowTumble) {
    for (int trial = 0; trial < 4; trial++) {
        Scenario s(3000 + trial, 1.0, 0.3, 0.1, 0.01);

        s.tumble(20000, 0.1 * kDegree);
        EXPECT_LT(s.errors().magnitude, 3e-3) << "trial " << trial;
    }
}

/* A noisy sensor */
TEST(Accuracy, HighNoise) {
    for (int trial = 0; trial < 10; trial++) {
        Scenario s(3100 + trial, 1.0, 0.3, 0.1, 0.05);

        s.tumble(3000);
        EXPECT_LT(s.errors().magnitude, 1.5e-2) << "trial " << trial;
    }
}

/*
With reference directions, a misalignment of the sensor is calibrated out as
well: the calibrated readings match the field as a vector.
*/
TEST(Accuracy, ReferenceCorrectsMisalignment) {
    for (int trial = 0; trial < 10; trial++) {
        Scenario s(3200 + trial, 1.0, 0.3, 0.1, 0.01);

        s.sensor.misalign(s.rng.direction(), 10.0 * kDegree);
        s.with_reference = true;
        s.tumble(3000);

        const Errors e = s.errors();
        EXPECT_LT(e.magnitude, 2.5e-3) << "trial " << trial;
        EXPECT_LT(e.vector, 5e-3) << "trial " << trial;
    }
}

/* Without a reference the misalignment is, necessarily, left in */
TEST(Accuracy, NoReferenceLeavesMisalignment) {
    Scenario s(3300, 1.0, 0.3, 0.1, 0.01);

    s.sensor.misalign(s.rng.direction(), 10.0 * kDegree);
    s.tumble(3000);

    const Errors e = s.errors();
    EXPECT_LT(e.magnitude, 2.5e-3);
    /*
    A rotation by 10 degrees moves a vector at right angles to its axis by
    2 sin(5 deg) = 0.174, and one along the axis not at all: over all
    directions that is an RMS of sqrt(2/3) * 0.174 = 0.142.
    */
    EXPECT_NEAR(0.142, e.vector, 0.01);
}

/* Errors in the reference average out */
TEST(Accuracy, ReferenceErrorsAverageOut) {
    for (int trial = 0; trial < 10; trial++) {
        Scenario s(3400 + trial, 1.0, 0.3, 0.1, 0.01);

        s.sensor.misalign(s.rng.direction(), 5.0 * kDegree);
        s.with_reference = true;
        s.reference_error = 3.0 * kDegree;
        s.tumble(3000);

        const Errors e = s.errors();
        EXPECT_LT(e.magnitude, 2.5e-3) << "trial " << trial;
        EXPECT_LT(e.vector, 8e-3) << "trial " << trial;
    }
}

/*
The bias and the shape of the scale estimate come from the measurements
alone: a reference, however bad, cannot affect them.
*/
TEST(Accuracy, ReferenceCannotAffectBiasOrShape) {
    Scenario with(3500, 1.0, 0.3, 0.1, 0.01);
    Scenario without(3500, 1.0, 0.3, 0.1, 0.01);

    with.with_reference = true;
    for (int i = 0; i < 2000; i++) {
        const Vec3 direction = with.motion.step(kDegree);
        const Vec3 reading = with.sensor.read(direction, with.rng);
        /* a reference that has nothing to do with the truth */
        const Vec3 nonsense = with.rng.direction();

        support::update(with.cal, reading, &nonsense);
        support::update(without.cal, reading);
    }

    const support::Calibration a = support::calibration_of(with.cal);
    const support::Calibration b = support::calibration_of(without.cal);
    const Mat3 ata = support::mul(support::transpose(a.a), a.a);
    const Mat3 btb = support::mul(support::transpose(b.a), b.a);

    for (int i = 0; i < 3; i++) {
        EXPECT_EQ(a.b[i], b.b[i]);
    }
    for (int i = 0; i < 9; i++) {
        EXPECT_NEAR(ata[i], btb[i], 1e-5);
    }
}

/*
The variances reported by TRICAL_estimate_get_ext are realistic: compare
them with the scatter of the estimates over many runs.
*/
TEST(Accuracy, ReportedVariances) {
    const int trials = 60;
    double sum_error2[TRICAL_STATE_DIM] = { 0 };
    double sum_variance[TRICAL_STATE_DIM] = { 0 };

    for (int trial = 0; trial < trials; trial++) {
        Scenario s(3600 + trial, 1.0, 0.0, 0.0, 0.01);
        float bias[3], scale[9], bias_variance[3], scale_variance[9];

        /* the same sensor every time, so the errors are comparable */
        s.sensor.b = Vec3{{ 0.2, -0.1, 0.3 }};
        s.sensor.a = Mat3{{ 1.1, 0.03, -0.02, 0.03, 0.95, 0.01, -0.02, 0.01,
                            1.05 }};
        s.scatter(1000);

        TRICAL_estimate_get_ext(&s.cal, bias, scale, bias_variance,
                                scale_variance);
        for (int i = 0; i < 3; i++) {
            const double error = bias[i] - s.sensor.b[i];

            sum_error2[i] += error * error;
            sum_variance[i] += bias_variance[i];
        }
        for (int i = 0; i < 9; i++) {
            const double error = scale[i] + ((i % 4 == 0) ? 1.0 : 0.0) -
                                 s.sensor.a[i];

            sum_error2[3 + i] += error * error;
            sum_variance[3 + i] += scale_variance[i];
        }
    }

    for (int i = 0; i < TRICAL_STATE_DIM; i++) {
        const double ratio = sum_variance[i] / sum_error2[i];

        EXPECT_GT(ratio, 0.5) << "state " << i;
        EXPECT_LT(ratio, 2.0) << "state " << i;
    }
}

/* Before the data determine the calibration, the variances say so */
TEST(Accuracy, ReportedVariancesReflectCoverage) {
    Scenario s(3700, 1.0, 0.3, 0.1, 0.01);
    float bias[3], scale[9], early_bias[3], early_scale[9], late_bias[3],
          late_scale[9];

    s.rest(500);
    TRICAL_estimate_get_ext(&s.cal, bias, scale, early_bias, early_scale);
    s.tumble(3000);
    TRICAL_estimate_get_ext(&s.cal, bias, scale, late_bias, late_scale);

    for (int i = 0; i < 3; i++) {
        EXPECT_GT(early_bias[i], 100 * late_bias[i]);
        EXPECT_LT(late_bias[i], 1e-5);
    }
    for (int i = 0; i < 9; i++) {
        EXPECT_GT(early_scale[i], 0.0f);
        EXPECT_LT(late_scale[i], 1e-5);
    }
}
