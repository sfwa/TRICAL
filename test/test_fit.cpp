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
#include <vector>

/* C++ has no 'restrict' qualifier. Just ignore it. */
#define restrict

#include "TRICAL.h"
#include "fit.h"
#include "support.h"

using support::Mat3;
using support::Random;
using support::Vec3;

namespace {

const bool kWideIsFloat = sizeof(trical_wide_t) == sizeof(float);
const double kEps = kWideIsFloat ? 3e-5 : 1e-11;

/* The quadric basis and its Jacobian, written out longhand. */
void basis(const double x[3], double xi[10]) {
    xi[0] = x[0] * x[0];
    xi[1] = x[1] * x[1];
    xi[2] = x[2] * x[2];
    xi[3] = 2 * x[1] * x[2];
    xi[4] = 2 * x[0] * x[2];
    xi[5] = 2 * x[0] * x[1];
    xi[6] = 2 * x[0];
    xi[7] = 2 * x[1];
    xi[8] = 2 * x[2];
    xi[9] = 1;
}

void jacobian(const double x[3], double j[10][3]) {
    const double rows[10][3] = {
        { 2 * x[0], 0, 0 }, { 0, 2 * x[1], 0 }, { 0, 0, 2 * x[2] },
        { 0, 2 * x[2], 2 * x[1] }, { 2 * x[2], 0, 2 * x[0] },
        { 2 * x[1], 2 * x[0], 0 },
        { 2, 0, 0 }, { 0, 2, 0 }, { 0, 0, 2 }, { 0, 0, 0 }
    };
    for (int i = 0; i < 10; i++) {
        for (int k = 0; k < 3; k++) {
            j[i][k] = rows[i][k];
        }
    }
}

std::vector<Vec3> random_points(Random &rng, int count, double spread) {
    std::vector<Vec3> points;
    for (int i = 0; i < count; i++) {
        points.push_back(Vec3{{ spread * rng.normal(), spread * rng.normal(),
                                spread * rng.normal() }});
    }
    return points;
}

void accumulate(const std::vector<Vec3> &points,
        trical_wide_t sums[TRICAL_NUM_MOMENTS], const Vec3 &offset) {
    for (int i = 0; i < TRICAL_NUM_MOMENTS; i++) {
        sums[i] = 0;
    }
    for (const Vec3 &p : points) {
        const float x[3] = { static_cast<float>(p[0] - offset[0]),
                             static_cast<float>(p[1] - offset[1]),
                             static_cast<float>(p[2] - offset[2]) };
        _trical_moments_accumulate(sums, x);
    }
}

/* Value of the quadric at x */
double evaluate(const trical_wide_t quadric[10], const Vec3 &x) {
    double xi[10], s = 0.0;
    basis(x.data(), xi);
    for (int i = 0; i < 10; i++) {
        s += xi[i] * static_cast<double>(quadric[i]);
    }
    return s;
}

} /* namespace */

/* The running sums are the sums of all monomials of degree <= 4 */
TEST(Fit, MomentsAccumulate) {
    Random rng(21);
    const std::vector<Vec3> points = random_points(rng, 40, 1.0);
    trical_wide_t sums[TRICAL_NUM_MOMENTS];
    int index = 0;

    accumulate(points, sums, Vec3{{ 0, 0, 0 }});

    /* graded lexicographic order: 1, x, y, z, xx, xy, xz, yy, ... */
    for (int degree = 0; degree <= 4; degree++) {
        for (int a = degree; a >= 0; a--) {
            for (int b = degree - a; b >= 0; b--) {
                const int c = degree - a - b;
                double expected = 0.0, magnitude = 0.0;

                for (const Vec3 &p : points) {
                    const float x = static_cast<float>(p[0]),
                                y = static_cast<float>(p[1]),
                                z = static_cast<float>(p[2]);
                    const double term = std::pow(x, a) * std::pow(y, b) *
                                        std::pow(z, c);
                    expected += term;
                    magnitude += std::fabs(term);
                }
                /* products are formed in single precision */
                EXPECT_NEAR(expected, sums[index], 1e-6 * magnitude + 1e-12)
                    << "x^" << a << " y^" << b << " z^" << c;
                index++;
            }
        }
    }
    EXPECT_EQ(TRICAL_NUM_MOMENTS, index);
}

/* Shifting the sums is the same as accumulating shifted points */
TEST(Fit, MomentsShift) {
    Random rng(22);
    const std::vector<Vec3> points = random_points(rng, 60, 1.0);
    const Vec3 offset = {{ 0.25, -0.5, 0.75 }};
    const trical_wide_t delta[3] = { static_cast<trical_wide_t>(offset[0]),
                                     static_cast<trical_wide_t>(offset[1]),
                                     static_cast<trical_wide_t>(offset[2]) };
    trical_wide_t shifted[TRICAL_NUM_MOMENTS], direct[TRICAL_NUM_MOMENTS];

    accumulate(points, shifted, Vec3{{ 0, 0, 0 }});
    _trical_moments_shift(shifted, delta);
    accumulate(points, direct, offset);

    for (int i = 0; i < TRICAL_NUM_MOMENTS; i++) {
        EXPECT_NEAR(direct[i], shifted[i], kWideIsFloat ? 2e-2 : 2e-4)
            << "moment " << i;
    }
}

/* The moment matrix is the mean of xi xi' */
/* The matrix whose product with a vector is given by `product` */
template <typename Product>
void matrix_of(Product product, trical_wide_t m[100]) {
    for (int j = 0; j < 10; j++) {
        trical_wide_t unit[10] = { 0 }, column[10];

        unit[j] = 1;
        product(unit, column);
        for (int i = 0; i < 10; i++) {
            m[i * 10 + j] = column[i];
        }
    }
}

void moment_matrix(const trical_wide_t sums[TRICAL_NUM_MOMENTS],
        trical_wide_t m[100]) {
    matrix_of([sums](const trical_wide_t *theta, trical_wide_t *y) {
        _trical_moment_product(sums, theta, y);
    }, m);
}

void normalisation_matrix(const trical_wide_t sums[TRICAL_NUM_MOMENTS],
        trical_wide_t n[100], int hyper) {
    trical_wide_t xi[10];

    _trical_basis_mean(sums, xi);
    matrix_of([&xi, hyper](const trical_wide_t *theta, trical_wide_t *y) {
        _trical_normalisation_product(xi, theta, hyper, y);
    }, n);
}

TEST(Fit, MomentMatrix) {
    Random rng(23);
    const std::vector<Vec3> points = random_points(rng, 50, 0.8);
    trical_wide_t sums[TRICAL_NUM_MOMENTS], m[100];
    double expected[100] = { 0 };

    accumulate(points, sums, Vec3{{ 0, 0, 0 }});
    moment_matrix(sums, m);

    for (const Vec3 &p : points) {
        const double x[3] = { static_cast<float>(p[0]),
                              static_cast<float>(p[1]),
                              static_cast<float>(p[2]) };
        double xi[10];

        basis(x, xi);
        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 10; j++) {
                expected[i * 10 + j] += xi[i] * xi[j] / points.size();
            }
        }
    }
    for (int i = 0; i < 100; i++) {
        EXPECT_NEAR(expected[i], m[i], 1e-5) << "entry " << i;
    }
}

/* The normalisation matrices are the means of J J' and of xi e' + e xi' */
TEST(Fit, NormalisationMatrix) {
    Random rng(24);
    const std::vector<Vec3> points = random_points(rng, 50, 0.8);
    trical_wide_t sums[TRICAL_NUM_MOMENTS], taubin[100], hyper[100];
    double expected[100] = { 0 }, mean_xi[10] = { 0 };

    accumulate(points, sums, Vec3{{ 0, 0, 0 }});
    normalisation_matrix(sums, taubin, 0);
    normalisation_matrix(sums, hyper, 1);

    for (const Vec3 &p : points) {
        const double x[3] = { static_cast<float>(p[0]),
                              static_cast<float>(p[1]),
                              static_cast<float>(p[2]) };
        double j[10][3], xi[10];

        jacobian(x, j);
        basis(x, xi);
        for (int a = 0; a < 10; a++) {
            mean_xi[a] += xi[a] / points.size();
            for (int b = 0; b < 10; b++) {
                expected[a * 10 + b] += (j[a][0] * j[b][0] +
                    j[a][1] * j[b][1] + j[a][2] * j[b][2]) / points.size();
            }
        }
    }

    for (int a = 0; a < 10; a++) {
        for (int b = 0; b < 10; b++) {
            const double e_a = (a < 3) ? 1.0 : 0.0, e_b = (b < 3) ? 1.0 : 0.0;

            EXPECT_NEAR(expected[a * 10 + b], taubin[a * 10 + b], 1e-5)
                << a << ", " << b;
            EXPECT_NEAR(expected[a * 10 + b] + mean_xi[a] * e_b +
                        e_a * mean_xi[b], hyper[a * 10 + b], 1e-5)
                << a << ", " << b;
        }
    }
}

/*
The noise correction: for noisy points the moment matrix is, on average, the
noise-free one plus sigma^2 times the "hyper" normalisation matrix. This is
the fact the whole method rests on.
*/
TEST(Fit, NoiseBiasOfMomentMatrix) {
    Random rng(25);
    const double sigma = 0.05;
    const int count = 400000;
    trical_wide_t clean[TRICAL_NUM_MOMENTS] = { 0 };
    trical_wide_t noisy[TRICAL_NUM_MOMENTS] = { 0 };
    trical_wide_t m_clean[100], m_noisy[100], n[100];

    for (int i = 0; i < count; i++) {
        const Vec3 d = rng.direction();
        const float x[3] = { static_cast<float>(d[0]),
                             static_cast<float>(d[1]),
                             static_cast<float>(d[2]) };
        const float y[3] = { static_cast<float>(d[0] + sigma * rng.normal()),
                             static_cast<float>(d[1] + sigma * rng.normal()),
                             static_cast<float>(d[2] + sigma * rng.normal()) };

        _trical_moments_accumulate(clean, x);
        _trical_moments_accumulate(noisy, y);
    }

    moment_matrix(clean, m_clean);
    moment_matrix(noisy, m_noisy);
    normalisation_matrix(clean, n, 1);

    for (int i = 0; i < 100; i++) {
        /* sampling error of the means, plus O(sigma^4) terms */
        EXPECT_NEAR(m_clean[i] + sigma * sigma * n[i], m_noisy[i], 1.5e-3)
            << "entry " << i;
    }
}

/* Quadric <-> calibration conversions are inverses of each other */
TEST(Fit, QuadricRoundTrip) {
    Random rng(26);

    for (int trial = 0; trial < 100; trial++) {
        const support::Sensor sensor = support::Sensor::random(rng, 1.0,
                                                               0.0, 0.3, 0.0);
        const Vec3 c = {{ rng.normal(), rng.normal(), rng.normal() }};
        const Mat3 ata = support::mul(support::transpose(sensor.a), sensor.a);
        trical_wide_t ata_w[9], centre_w[3], quadric[10];
        trical_wide_t centre_out[3], a_out[9], a_inverse[9];

        for (int i = 0; i < 9; i++) {
            ata_w[i] = static_cast<trical_wide_t>(ata[i]);
        }
        for (int i = 0; i < 3; i++) {
            centre_w[i] = static_cast<trical_wide_t>(c[i]);
        }

        _trical_quadric_from_calibration(ata_w, centre_w, quadric);

        /* unit norm, and zero on the ellipsoid */
        double norm2 = 0.0;
        for (int i = 0; i < 10; i++) {
            norm2 += static_cast<double>(quadric[i]) * quadric[i];
        }
        EXPECT_NEAR(1.0, norm2, 10 * kEps);
        for (int k = 0; k < 5; k++) {
            const Vec3 on_surface = support::mul(support::inverse(sensor.a),
                                                 rng.direction());
            const Vec3 x = {{ on_surface[0] + c[0], on_surface[1] + c[1],
                              on_surface[2] + c[2] }};
            EXPECT_NEAR(0.0, evaluate(quadric, x), 30 * kEps);
        }

        ASSERT_EQ(0, _trical_quadric_to_calibration(quadric, centre_out,
                                                    a_out, a_inverse));
        for (int i = 0; i < 3; i++) {
            EXPECT_NEAR(c[i], centre_out[i], 100 * kEps);
        }
        /* sensor.a is symmetric positive definite, so it is recovered */
        const Mat3 identity = support::identity();
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                double product = 0.0;
                for (int k = 0; k < 3; k++) {
                    product += static_cast<double>(a_out[i * 3 + k]) *
                               a_inverse[k * 3 + j];
                }
                EXPECT_NEAR(sensor.a[i * 3 + j], a_out[i * 3 + j],
                            100 * kEps);
                EXPECT_NEAR(identity[i * 3 + j], product, 100 * kEps);
            }
        }
    }
}

/* The sign and scale of a quadric do not matter */
TEST(Fit, QuadricToCalibrationIgnoresScale) {
    const trical_wide_t ata[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
    const trical_wide_t centre[3] = { 0.5, -0.25, 0.125 };
    trical_wide_t quadric[10], scaled[10], c1[3], c2[3], a1[9], a2[9];

    _trical_quadric_from_calibration(ata, centre, quadric);
    for (int i = 0; i < 10; i++) {
        scaled[i] = static_cast<trical_wide_t>(-7.5) * quadric[i];
    }

    ASSERT_EQ(0, _trical_quadric_to_calibration(quadric, c1, a1, nullptr));
    ASSERT_EQ(0, _trical_quadric_to_calibration(scaled, c2, a2, nullptr));
    for (int i = 0; i < 3; i++) {
        EXPECT_NEAR(c1[i], c2[i], 10 * kEps);
    }
    for (int i = 0; i < 9; i++) {
        EXPECT_NEAR(a1[i], a2[i], 10 * kEps);
    }
}

/* Quadrics that are not credible ellipsoids are rejected */
TEST(Fit, QuadricToCalibrationRejectsNonEllipsoids) {
    /* (G11, G22, G33, G23, G13, G12, q1, q2, q3, q0) */
    const trical_wide_t hyperboloid[10] = { 1, 1, -1, 0, 0, 0, 0, 0, 0, -1 };
    const trical_wide_t cylinder[10] = { 1, 1, 0, 0, 0, 0, 0, 0, 0, -1 };
    const trical_wide_t plane[10] = { 0, 0, 0, 0, 0, 0, 1, 0, 0, -1 };
    const trical_wide_t empty[10] = { 1, 1, 1, 0, 0, 0, 0, 0, 0, 1 };
    const trical_wide_t needle[10] = { 1, 1, 0.01f, 0, 0, 0, 0, 0, 0, -1 };
    const trical_wide_t huge[10] = { 1, 1, 1, 0, 0, 0, 0, 0, 0, -100 };
    const trical_wide_t tiny[10] = { 100, 100, 100, 0, 0, 0, 0, 0, 0, -1 };
    const trical_wide_t zero[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    const trical_wide_t sphere[10] = { 1, 1, 1, 0, 0, 0, 0, 0, 0, -1 };
    trical_wide_t centre[3], a[9];

    EXPECT_EQ(-1, _trical_quadric_to_calibration(hyperboloid, centre, a,
                                                 nullptr));
    EXPECT_EQ(-1, _trical_quadric_to_calibration(cylinder, centre, a,
                                                 nullptr));
    EXPECT_EQ(-1, _trical_quadric_to_calibration(plane, centre, a, nullptr));
    EXPECT_EQ(-1, _trical_quadric_to_calibration(empty, centre, a, nullptr));
    EXPECT_EQ(-1, _trical_quadric_to_calibration(needle, centre, a, nullptr));
    EXPECT_EQ(-1, _trical_quadric_to_calibration(huge, centre, a, nullptr));
    EXPECT_EQ(-1, _trical_quadric_to_calibration(tiny, centre, a, nullptr));
    EXPECT_EQ(-1, _trical_quadric_to_calibration(zero, centre, a, nullptr));
    EXPECT_EQ(0, _trical_quadric_to_calibration(sphere, centre, a, nullptr));
}

/* Translating a quadric moves its zero set */
TEST(Fit, QuadricTranslate) {
    Random rng(27);
    const trical_wide_t ata[9] = { 2, 0.3f, 0, 0.3f, 1, -0.2f, 0, -0.2f,
                                   1.5f };
    const trical_wide_t centre[3] = { 0.3f, 0.2f, -0.1f };
    const trical_wide_t delta[3] = { 1.5f, -2.0f, 0.5f };
    trical_wide_t quadric[10], moved[10];

    _trical_quadric_from_calibration(ata, centre, quadric);
    for (int i = 0; i < 10; i++) {
        moved[i] = quadric[i];
    }
    _trical_quadric_translate(moved, delta);

    /* the same ellipsoid, with its centre at centre - delta */
    trical_wide_t c[3], a[9], c0[3], a0[9];
    ASSERT_EQ(0, _trical_quadric_to_calibration(quadric, c0, a0, nullptr));
    ASSERT_EQ(0, _trical_quadric_to_calibration(moved, c, a, nullptr));
    for (int i = 0; i < 3; i++) {
        EXPECT_NEAR(c0[i] - delta[i], c[i], 100 * kEps);
    }
    for (int i = 0; i < 9; i++) {
        EXPECT_NEAR(a0[i], a[i], 100 * kEps);
    }
}

/* The polar rotation of (rotation x symmetric positive definite) */
TEST(Fit, PolarRotation) {
    Random rng(28);

    for (int trial = 0; trial < 100; trial++) {
        const Mat3 r = support::rotation(rng.direction(),
                                         rng.uniform(-3.0, 3.0));
        const support::Sensor shape = support::Sensor::random(rng, 1.0, 0.0,
                                                              0.5, 0.0);
        const Mat3 k = support::mul(r, shape.a);
        trical_wide_t k_w[9], q[9];

        for (int i = 0; i < 9; i++) {
            k_w[i] = static_cast<trical_wide_t>(k[i]);
        }
        ASSERT_EQ(0, _trical_polar_rotation(k_w, q));
        for (int i = 0; i < 9; i++) {
            EXPECT_NEAR(r[i], q[i], 100 * kEps);
        }
    }
}

/* The result is a proper rotation even if the input is a reflection */
TEST(Fit, PolarRotationOfReflection) {
    const trical_wide_t k[9] = { 1, 0, 0, 0, 1, 0, 0, 0, -0.1f };
    trical_wide_t q[9];

    ASSERT_EQ(0, _trical_polar_rotation(k, q));
    EXPECT_NEAR(1.0, trical_det3(q), 100 * kEps);
    /* the best rotation leaves the well-determined axes alone */
    EXPECT_NEAR(1.0, q[0], 100 * kEps);
    EXPECT_NEAR(1.0, q[4], 100 * kEps);
    EXPECT_NEAR(1.0, q[8], 100 * kEps);
}

/*
Two independent directions are enough to fix a rotation: reference
directions that all lie in a plane, say. The result is a proper rotation to
within rounding however small the third singular value is.
*/
TEST(Fit, PolarRotationOfRankTwo) {
    Random rng(31);

    for (int trial = 0; trial < 200; trial++) {
        const Mat3 u = support::rotation(rng.direction(),
                                         rng.uniform(-3.0, 3.0));
        const Mat3 v = support::rotation(rng.direction(),
                                         rng.uniform(-3.0, 3.0));
        const double third[4] = { 0.0, 1e-3, 1e-7, -1e-3 };
        const double s[3] = { rng.uniform(0.3, 1.0), rng.uniform(0.3, 1.0),
                              third[trial % 4] };
        trical_wide_t k[9], q[9];

        /* k = u diag(s) v' */
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                double sum = 0.0;
                for (int l = 0; l < 3; l++) {
                    sum += u[i * 3 + l] * s[l] * v[j * 3 + l];
                }
                k[i * 3 + j] = static_cast<trical_wide_t>(sum);
            }
        }

        ASSERT_EQ(0, _trical_polar_rotation(k, q)) << "trial " << trial;

        /* orthogonal, with determinant +1 */
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                double dot = 0.0;
                for (int l = 0; l < 3; l++) {
                    dot += static_cast<double>(q[l * 3 + i]) * q[l * 3 + j];
                }
                EXPECT_NEAR((i == j) ? 1.0 : 0.0, dot, 100 * kEps)
                    << "trial " << trial;
            }
        }
        EXPECT_NEAR(1.0, trical_det3(q), 100 * kEps) << "trial " << trial;

        /*
        and it is u v', whatever the sign of the third singular value (u and
        v being rotations)
        */
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                double expected = 0.0;
                for (int l = 0; l < 3; l++) {
                    expected += u[i * 3 + l] * v[j * 3 + l];
                }
                EXPECT_NEAR(expected, q[i * 3 + j], 2000 * kEps)
                    << "trial " << trial;
            }
        }
    }
}

/* Fewer than two independent directions are not, and that is reported */
TEST(Fit, PolarRotationRejectsRankDeficient) {
    const trical_wide_t rank1[9] = { 1, 2, 3, 2, 4, 6, 3, 6, 9 };
    const trical_wide_t zero[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    trical_wide_t q[9];

    EXPECT_EQ(-1, _trical_polar_rotation(rank1, q));
    EXPECT_EQ(-1, _trical_polar_rotation(zero, q));
}

/*
Exact data on an ellipsoid are fitted exactly, wherever the ellipsoid is and
whatever the units.
*/
TEST(Fit, RefitRecoversEllipsoidExactly) {
    Random rng(29);
    const double norms[3] = { 1.0, 47.0, 0.003 };

    for (int trial = 0; trial < 60; trial++) {
        const double field_norm = norms[trial % 3];
        const support::Sensor sensor = support::Sensor::random(rng,
            field_norm, rng.uniform(0.0, 4.0), 0.3, 0.0);
        TRICAL_instance_t cal;

        TRICAL_init(&cal);
        TRICAL_norm_set(&cal, static_cast<float>(field_norm));
        for (int i = 0; i < 400; i++) {
            support::update(cal, sensor.read(rng.direction(), rng));
        }

        const support::Errors errors = support::errors_of(cal, sensor);
        EXPECT_LT(errors.magnitude, kWideIsFloat ? 2e-3 : 2e-5)
            << "trial " << trial;
        EXPECT_LT(errors.vector, kWideIsFloat ? 2e-3 : 2e-5)
            << "trial " << trial;
    }
}

/* The result does not depend on the order of the measurements */
TEST(Fit, OrderIndependent) {
    Random rng(30);
    const support::Sensor sensor = support::Sensor::random(rng, 1.0, 0.5,
                                                           0.15, 0.01);
    std::vector<Vec3> readings;
    TRICAL_instance_t forward, backward;

    /*
    Well separated, so that every reading starts a new cell, and none is in
    the cell of the reading two before or after it -- which would make the
    one in between look like a glitch, in one order or the other. Each end
    of the list is an exception: that is where the list is closed, below.
    */
    const auto separate = [](const Vec3 &p, const Vec3 &q) {
        double distance2 = 0;

        for (int k = 0; k < 3; k++) {
            distance2 += (p[k] - q[k]) * (p[k] - q[k]);
        }
        return distance2 > 0.5 * 0.5;
    };

    while (readings.size() < 300) {
        const Vec3 reading = sensor.read(rng.direction(), rng);
        const size_t n = readings.size();

        if ((n < 1 || separate(reading, readings[n - 1])) &&
                (n < 2 || separate(reading, readings[n - 2]))) {
            readings.push_back(reading);
        }
    }

    /*
    A reading is used once the next one has arrived, so finish each pass
    with a reading that is not part of the comparison: one that is well
    away from the readings at both ends of the list.
    */
    Vec3 last;
    do {
        last = sensor.read(rng.direction(), rng);
    } while (!separate(last, readings[0]) || !separate(last, readings[1]) ||
             !separate(last, readings[298]) ||
             !separate(last, readings[299]));

    TRICAL_init(&forward);
    TRICAL_init(&backward);
    for (size_t i = 0; i < readings.size(); i++) {
        support::update(forward, readings[i]);
        support::update(backward, readings[readings.size() - 1 - i]);
    }
    support::update(forward, last);
    support::update(backward, last);
    ASSERT_EQ(300.0, static_cast<double>(forward.priv.moments[0]));
    ASSERT_EQ(300.0, static_cast<double>(backward.priv.moments[0]));

    /* the sums are the same whatever the order, up to rounding */
    for (int i = 0; i < TRICAL_NUM_MOMENTS; i++) {
        EXPECT_NEAR(forward.priv.moments[i], backward.priv.moments[i],
                    kWideIsFloat ? 2e-2 : 1e-9) << "sum " << i;
    }

    /* bring both estimates up to date with their sums */
    _trical_refit(&forward);
    _trical_refit(&backward);

    /*
    The noise level that sets the weight of the prior is refined from one
    refit to the next, though, so the estimates agree to a small fraction of
    their statistical error (about 1e-3 here) rather than exactly.
    */
    const support::Calibration f = support::calibration_of(forward);
    const support::Calibration b = support::calibration_of(backward);
    for (int i = 0; i < 3; i++) {
        EXPECT_NEAR(f.b[i], b.b[i], 1e-4);
    }
    for (int i = 0; i < 9; i++) {
        EXPECT_NEAR(f.a[i], b.a[i], 1e-4);
    }
}
