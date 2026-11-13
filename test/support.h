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

#ifndef TRICAL_TEST_SUPPORT_H
#define TRICAL_TEST_SUPPORT_H

/*
Support code for the tests: a simulated sensor, error metrics, and an
independent batch maximum-likelihood fit to measure TRICAL against.
Everything here is double precision and deliberately straightforward.
*/

#include <array>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

#include "TRICAL.h"

namespace support {

typedef std::array<double, 3> Vec3;
typedef std::array<double, 9> Mat3; /* row-major */

/*
Random numbers that are the same on every platform: the standard fixes the
output of mt19937 but not of the distributions, so those are done by hand.
*/
class Random {
public:
    explicit Random(uint32_t seed) : gen_(seed), have_spare_(false),
        spare_(0.0) {}

    /* uniform on (0, 1) */
    double uniform() {
        return (static_cast<double>(gen_()) + 0.5) / 4294967296.0;
    }

    double uniform(double lo, double hi) {
        return lo + (hi - lo) * uniform();
    }

    /* standard normal (Box-Muller) */
    double normal() {
        if (have_spare_) {
            have_spare_ = false;
            return spare_;
        }
        const double r = std::sqrt(-2.0 * std::log(uniform()));
        const double t = 6.283185307179586 * uniform();
        spare_ = r * std::sin(t);
        have_spare_ = true;
        return r * std::cos(t);
    }

    /* uniformly distributed unit vector */
    Vec3 direction() {
        for (;;) {
            Vec3 v = {{ normal(), normal(), normal() }};
            const double n = std::sqrt(v[0] * v[0] + v[1] * v[1] +
                                       v[2] * v[2]);
            if (n > 1e-6) {
                return Vec3{{ v[0] / n, v[1] / n, v[2] / n }};
            }
        }
    }

private:
    std::mt19937 gen_;
    bool have_spare_;
    double spare_;
};

inline Mat3 identity() {
    return Mat3{{ 1, 0, 0, 0, 1, 0, 0, 0, 1 }};
}

inline Mat3 mul(const Mat3 &a, const Mat3 &b) {
    Mat3 c;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            c[i * 3 + j] = a[i * 3] * b[j] + a[i * 3 + 1] * b[3 + j] +
                           a[i * 3 + 2] * b[6 + j];
        }
    }
    return c;
}

inline Vec3 mul(const Mat3 &a, const Vec3 &v) {
    return Vec3{{ a[0] * v[0] + a[1] * v[1] + a[2] * v[2],
                  a[3] * v[0] + a[4] * v[1] + a[5] * v[2],
                  a[6] * v[0] + a[7] * v[1] + a[8] * v[2] }};
}

inline Mat3 transpose(const Mat3 &a) {
    return Mat3{{ a[0], a[3], a[6], a[1], a[4], a[7], a[2], a[5], a[8] }};
}

inline double det(const Mat3 &a) {
    return a[0] * (a[4] * a[8] - a[5] * a[7]) -
           a[1] * (a[3] * a[8] - a[5] * a[6]) +
           a[2] * (a[3] * a[7] - a[4] * a[6]);
}

inline Mat3 inverse(const Mat3 &a) {
    const double d = 1.0 / det(a);
    return Mat3{{
        (a[4] * a[8] - a[5] * a[7]) * d, (a[2] * a[7] - a[1] * a[8]) * d,
        (a[1] * a[5] - a[2] * a[4]) * d,
        (a[5] * a[6] - a[3] * a[8]) * d, (a[0] * a[8] - a[2] * a[6]) * d,
        (a[2] * a[3] - a[0] * a[5]) * d,
        (a[3] * a[7] - a[4] * a[6]) * d, (a[1] * a[6] - a[0] * a[7]) * d,
        (a[0] * a[4] - a[1] * a[3]) * d }};
}

inline double norm(const Vec3 &v) {
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

/* Rotation by `angle` radians about the unit vector `axis` (Rodrigues). */
inline Mat3 rotation(const Vec3 &axis, double angle) {
    const double c = std::cos(angle), s = std::sin(angle), t = 1.0 - c;
    const double x = axis[0], y = axis[1], z = axis[2];
    return Mat3{{ t * x * x + c, t * x * y - s * z, t * x * z + s * y,
                  t * x * y + s * z, t * y * y + c, t * y * z - s * x,
                  t * x * z - s * y, t * y * z + s * x, t * z * z + c }};
}

/*
A simulated sensor: a field H in the sensor frame reads as
    B = a^-1 H + b + noise
so the correct calibration is H = a (B - b). `a` is symmetric positive
definite unless `misalign` is called.
*/
struct Sensor {
    Mat3 a;
    Vec3 b;
    double field_norm;
    double noise;

    Sensor() : a(identity()), b{{ 0, 0, 0 }}, field_norm(1.0), noise(0.0) {}

    /*
    Random bias of magnitude `bias` field norms, and random symmetric scale
    matrix with principal scale factors within 1 +/- `scale_error`.
    */
    static Sensor random(Random &rng, double field_norm, double bias,
            double scale_error, double noise) {
        Sensor s;
        const Vec3 direction = rng.direction();

        /* (one at a time: the order of evaluation matters) */
        const double angle = rng.uniform(0.0, 3.141592653589793);
        const Vec3 axis = rng.direction();
        const Mat3 v = rotation(axis, angle);
        Mat3 d = identity();

        for (int i = 0; i < 3; i++) {
            d[i * 3 + i] = 1.0 + scale_error * rng.uniform(-1.0, 1.0);
            s.b[i] = direction[i] * bias * field_norm;
        }
        s.a = mul(mul(v, d), transpose(v));
        s.field_norm = field_norm;
        s.noise = noise * field_norm;
        return s;
    }

    /* Rotates the sensor axes relative to the frame of the field. */
    void misalign(const Vec3 &axis, double angle) {
        a = mul(rotation(axis, angle), a);
    }

    /* Reading for a field of magnitude field_norm along `direction`. */
    Vec3 read(const Vec3 &direction, Random &rng) const {
        const Vec3 h = {{ direction[0] * field_norm,
                          direction[1] * field_norm,
                          direction[2] * field_norm }};
        Vec3 m = mul(inverse(a), h);
        for (int i = 0; i < 3; i++) {
            m[i] += b[i] + noise * rng.normal();
        }
        return m;
    }
};

/*
Field direction in the frame of a sensor that tumbles smoothly: each step
rotates by `step` radians about an axis that wanders slowly (or about a
fixed axis, which sweeps out a cone -- "planar" motion).
*/
class Tumble {
public:
    explicit Tumble(Random &rng) : rng_(rng), direction_(rng.direction()),
        axis_(rng.direction()) {}

    const Vec3 &direction() const {
        return direction_;
    }

    const Vec3 &step(double step) {
        for (int i = 0; i < 3; i++) {
            axis_[i] += 0.05 * rng_.normal();
        }
        normalise(axis_);
        return rotate(axis_, step);
    }

    const Vec3 &step_about(const Vec3 &axis, double step) {
        Vec3 a = axis;
        normalise(a);
        return rotate(a, step);
    }

private:
    static void normalise(Vec3 &v) {
        const double n = norm(v);
        for (int i = 0; i < 3; i++) {
            v[i] /= n;
        }
    }

    const Vec3 &rotate(const Vec3 &axis, double angle) {
        direction_ = mul(rotation(axis, angle), direction_);
        normalise(direction_);
        return direction_;
    }

    Random &rng_;
    Vec3 direction_;
    Vec3 axis_;
};

/* A calibration in TRICAL's form: calibrated = (I + d) (m - b). */
struct Calibration {
    Vec3 b;
    Mat3 a; /* I + D */
};

inline Calibration calibration_of(const TRICAL_instance_t &instance) {
    float bias[3], scale[9];
    Calibration c;

    TRICAL_estimate_get(&instance, bias, scale);
    for (int i = 0; i < 3; i++) {
        c.b[i] = bias[i];
    }
    for (int i = 0; i < 9; i++) {
        c.a[i] = scale[i] + ((i % 4 == 0) ? 1.0 : 0.0);
    }
    return c;
}

inline bool is_finite(const Calibration &c) {
    for (int i = 0; i < 3; i++) {
        if (!std::isfinite(c.b[i])) {
            return false;
        }
    }
    for (int i = 0; i < 9; i++) {
        if (!std::isfinite(c.a[i])) {
            return false;
        }
    }
    return true;
}

/*
Errors of a calibration over the whole sphere of field directions (a
Fibonacci lattice), as fractions of the field norm:
- magnitude: RMS of |calibrated| - field_norm
- vector: RMS of |calibrated - field|
*/
struct Errors {
    double magnitude;
    double vector;
};

inline Errors errors_of(const Calibration &c, const Sensor &sensor,
        int points = 500) {
    const Mat3 forward = inverse(sensor.a);
    double sum_m = 0.0, sum_v = 0.0;

    for (int k = 0; k < points; k++) {
        const double z = 1.0 - (2.0 * k + 1.0) / points;
        const double r = std::sqrt(1.0 - z * z);
        const double phi = 2.399963229728653 * k; /* golden angle */
        const Vec3 h = {{ sensor.field_norm * r * std::cos(phi),
                          sensor.field_norm * r * std::sin(phi),
                          sensor.field_norm * z }};
        Vec3 m = mul(forward, h);
        for (int i = 0; i < 3; i++) {
            m[i] += sensor.b[i] - c.b[i];
        }
        const Vec3 out = mul(c.a, m);
        const double magnitude = norm(out) - sensor.field_norm;
        const Vec3 diff = {{ out[0] - h[0], out[1] - h[1], out[2] - h[2] }};

        sum_m += magnitude * magnitude;
        sum_v += diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
    }

    Errors e;
    e.magnitude = std::sqrt(sum_m / points) / sensor.field_norm;
    e.vector = std::sqrt(sum_v / points) / sensor.field_norm;
    return e;
}

inline Errors errors_of(const TRICAL_instance_t &instance,
        const Sensor &sensor) {
    return errors_of(calibration_of(instance), sensor);
}

/* Feeds a reading to TRICAL, with an optional reference direction. */
inline void update(TRICAL_instance_t &instance, const Vec3 &m,
        const Vec3 *reference = nullptr) {
    float measurement[3] = { static_cast<float>(m[0]),
                             static_cast<float>(m[1]),
                             static_cast<float>(m[2]) };

    if (reference) {
        float field[3] = { static_cast<float>((*reference)[0]),
                           static_cast<float>((*reference)[1]),
                           static_cast<float>((*reference)[2]) };
        TRICAL_estimate_update(&instance, measurement, field);
    } else {
        TRICAL_estimate_update(&instance, measurement, nullptr);
    }
}

/*
Batch maximum-likelihood fit of a bias and a symmetric scale matrix to
`samples`, by Levenberg-Marquardt on the distance of each sample from the
ellipsoid |a (m - b)| = field_norm (to first order: the Sampson distance).
Started from `start`; this is the yardstick TRICAL is measured against, not
something a small processor would want to run.
*/
namespace detail {

inline void unpack(const double p[9], Calibration &c) {
    c.b = Vec3{{ p[0], p[1], p[2] }};
    c.a = Mat3{{ p[3], p[4], p[5], p[4], p[6], p[7], p[5], p[7], p[8] }};
}

inline double distance(const Calibration &c, const Vec3 &m,
        double field_norm) {
    const Vec3 v = {{ m[0] - c.b[0], m[1] - c.b[1], m[2] - c.b[2] }};
    const Vec3 out = mul(c.a, v);
    const double n = norm(out);
    const Vec3 unit = {{ out[0] / n, out[1] / n, out[2] / n }};
    /* gradient of |a (m - b)| with respect to m is a' unit */
    const Vec3 gradient = mul(transpose(c.a), unit);

    return (n - field_norm) / norm(gradient);
}

inline double cost(const double p[9], const std::vector<Vec3> &samples,
        double field_norm) {
    Calibration c;
    double sum = 0.0;

    unpack(p, c);
    for (const Vec3 &m : samples) {
        const double d = distance(c, m, field_norm);
        sum += d * d;
    }
    return sum;
}

/* Solves the symmetric positive definite 9x9 system a x = b in place. */
inline bool solve9(double a[81], double b[9]) {
    for (int j = 0; j < 9; j++) {
        for (int k = 0; k < j; k++) {
            a[j * 9 + j] -= a[j * 9 + k] * a[j * 9 + k];
        }
        if (!(a[j * 9 + j] > 0.0)) {
            return false;
        }
        a[j * 9 + j] = std::sqrt(a[j * 9 + j]);
        for (int i = j + 1; i < 9; i++) {
            for (int k = 0; k < j; k++) {
                a[i * 9 + j] -= a[i * 9 + k] * a[j * 9 + k];
            }
            a[i * 9 + j] /= a[j * 9 + j];
        }
    }
    for (int i = 0; i < 9; i++) {
        for (int k = 0; k < i; k++) {
            b[i] -= a[i * 9 + k] * b[k];
        }
        b[i] /= a[i * 9 + i];
    }
    for (int i = 8; i >= 0; i--) {
        for (int k = i + 1; k < 9; k++) {
            b[i] -= a[k * 9 + i] * b[k];
        }
        b[i] /= a[i * 9 + i];
    }
    return true;
}

} /* namespace detail */

inline Calibration maximum_likelihood(const std::vector<Vec3> &samples,
        double field_norm, const Calibration &start) {
    double p[9] = { start.b[0], start.b[1], start.b[2],
                    start.a[0], start.a[1], start.a[2],
                    start.a[4], start.a[5], start.a[8] };
    double best = detail::cost(p, samples, field_norm);
    double damping = 1e-6;

    for (int iteration = 0; iteration < 60; iteration++) {
        double jtj[81] = { 0 }, jtr[9] = { 0 };

        for (const Vec3 &m : samples) {
            Calibration c;
            double row[9];

            detail::unpack(p, c);
            const double r = detail::distance(c, m, field_norm);
            for (int k = 0; k < 9; k++) {
                const double h = 1e-6 * (k < 3 ? field_norm : 1.0);
                double q[9];
                Calibration plus, minus;

                std::copy(p, p + 9, q);
                q[k] = p[k] + h;
                detail::unpack(q, plus);
                q[k] = p[k] - h;
                detail::unpack(q, minus);
                row[k] = (detail::distance(plus, m, field_norm) -
                          detail::distance(minus, m, field_norm)) / (2 * h);
            }
            for (int i = 0; i < 9; i++) {
                jtr[i] += row[i] * r;
                for (int j = 0; j < 9; j++) {
                    jtj[i * 9 + j] += row[i] * row[j];
                }
            }
        }

        bool improved = false;
        for (int attempt = 0; attempt < 10 && !improved; attempt++) {
            double a[81], step[9], q[9];

            std::copy(jtj, jtj + 81, a);
            for (int i = 0; i < 9; i++) {
                a[i * 9 + i] *= 1.0 + damping;
                step[i] = -jtr[i];
            }
            if (detail::solve9(a, step)) {
                for (int i = 0; i < 9; i++) {
                    q[i] = p[i] + step[i];
                }
                const double trial = detail::cost(q, samples, field_norm);
                if (std::isfinite(trial) && trial <= best) {
                    const double gain = best - trial;

                    std::copy(q, q + 9, p);
                    best = trial;
                    damping = std::max(damping * 0.1, 1e-12);
                    improved = true;
                    if (gain <= 1e-14 * best) {
                        iteration = 1000; /* converged */
                    }
                    break;
                }
            }
            damping *= 10.0;
        }
        if (!improved) {
            break;
        }
    }

    Calibration result;
    detail::unpack(p, result);
    return result;
}

/* The correct calibration of `sensor`, as a starting point for the above. */
inline Calibration truth_of(const Sensor &sensor) {
    Calibration c;
    c.b = sensor.b;
    c.a = sensor.a;
    return c;
}

} /* namespace support */

#endif
