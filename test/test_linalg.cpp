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

#include <algorithm>
#include <cmath>
#include <limits>

/* C++ has no 'restrict' qualifier. Just ignore it. */
#define restrict

#include "TRICAL.h"
#include "linalg.h"
#include "support.h"

namespace {

/* Tolerance appropriate to the wide type under test. */
const double kEps = (sizeof(trical_wide_t) == sizeof(float)) ? 2e-5 : 1e-12;

} /* namespace */

/*
Test LLT decomposition of a square positive definite matrix via the Cholesky
method. The factor is returned in the lower triangle.
*/
TEST(LinAlg, CholeskyLLT) {
    trical_wide_t m[16] = {
        18, 22,  54,  42,
        22, 70,  86,  62,
        54, 86, 174, 134,
        42, 62, 134, 106
    };
    const double expected[16] = {
        4.24264, 0, 0, 0,
        5.18545, 6.5659051, 0, 0,
        12.72792, 3.0460384, 1.6497422, 0,
        9.8994951, 1.6245539, 1.8497111, 1.3926213
    };

    ASSERT_EQ(0, trical_cholesky(m, 4));
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j <= i; j++) {
            EXPECT_NEAR(expected[i * 4 + j], m[i * 4 + j], 2e-5);
        }
    }
}

/* Only the lower triangle of the input is read */
TEST(LinAlg, CholeskyReadsLowerTriangleOnly) {
    trical_wide_t m[9] = {
        4, 99, 99,
        2, 5, 99,
        2, 3, 6
    };

    ASSERT_EQ(0, trical_cholesky(m, 3));
    EXPECT_NEAR(2.0, m[0], kEps);
    EXPECT_NEAR(1.0, m[3], kEps);
    EXPECT_NEAR(2.0, m[4], kEps);
    EXPECT_NEAR(1.0, m[6], kEps);
    EXPECT_NEAR(1.0, m[7], kEps);
    EXPECT_NEAR(2.0, m[8], kEps);
}

/* Matrices that are not positive definite are reported, not factorised */
TEST(LinAlg, CholeskyRejectsIndefinite) {
    trical_wide_t indefinite[4] = { 1, 2, 2, 1 };
    trical_wide_t singular[4] = { 1, 1, 1, 1 };
    trical_wide_t zero[4] = { 0, 0, 0, 0 };
    trical_wide_t nan[4] = { 1, 0, 0,
                             std::numeric_limits<trical_wide_t>::quiet_NaN() };

    EXPECT_EQ(-1, trical_cholesky(indefinite, 2));
    EXPECT_EQ(-1, trical_cholesky(singular, 2));
    EXPECT_EQ(-1, trical_cholesky(zero, 2));
    EXPECT_EQ(-1, trical_cholesky(nan, 2));
}

/* Solve a random symmetric positive definite 10x10 system */
TEST(LinAlg, CholeskySolve) {
    support::Random rng(11);
    trical_wide_t a[100], l[100], x[10], b[10];

    /* a = g g' + I */
    double g[100];
    for (int i = 0; i < 100; i++) {
        g[i] = rng.normal();
    }
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            double s = (i == j) ? 1.0 : 0.0;
            for (int k = 0; k < 10; k++) {
                s += g[i * 10 + k] * g[j * 10 + k];
            }
            a[i * 10 + j] = static_cast<trical_wide_t>(s);
            l[i * 10 + j] = a[i * 10 + j];
        }
    }
    for (int i = 0; i < 10; i++) {
        b[i] = static_cast<trical_wide_t>(rng.normal());
        x[i] = b[i];
    }

    ASSERT_EQ(0, trical_cholesky(l, 10));
    trical_cholesky_solve(l, 10, x);

    for (int i = 0; i < 10; i++) {
        double s = 0.0;
        for (int j = 0; j < 10; j++) {
            s += static_cast<double>(a[i * 10 + j]) * x[j];
        }
        EXPECT_NEAR(b[i], s, 200 * kEps);
    }
}

/* Eigen-decomposition of symmetric 3x3 matrices: a = v diag(w) v' */
TEST(LinAlg, EigenSym3) {
    support::Random rng(12);

    for (int trial = 0; trial < 200; trial++) {
        trical_wide_t a[9], w[3], v[9];
        double scale = std::pow(10.0, rng.uniform(-3.0, 3.0));

        for (int i = 0; i < 3; i++) {
            for (int j = i; j < 3; j++) {
                a[i * 3 + j] = static_cast<trical_wide_t>(scale *
                                                          rng.normal());
                a[j * 3 + i] = a[i * 3 + j];
            }
        }

        trical_eigen_sym3(a, w, v);

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                double vtv = 0.0, vwv = 0.0;
                for (int k = 0; k < 3; k++) {
                    vtv += static_cast<double>(v[k * 3 + i]) * v[k * 3 + j];
                    vwv += static_cast<double>(v[i * 3 + k]) * w[k] *
                           v[j * 3 + k];
                }
                EXPECT_NEAR((i == j) ? 1.0 : 0.0, vtv, 50 * kEps);
                EXPECT_NEAR(a[i * 3 + j], vwv, 50 * kEps * scale);
            }
        }
    }
}

/* Only the upper triangle of the input is read, and it is not modified */
TEST(LinAlg, EigenSym3ReadsUpperTriangleOnly) {
    const trical_wide_t a[9] = {
        2, 1, 0,
        99, 2, 0,
        99, 99, 5
    };
    trical_wide_t w[3], v[9];

    trical_eigen_sym3(a, w, v);

    double sorted[3] = { static_cast<double>(w[0]),
                         static_cast<double>(w[1]),
                         static_cast<double>(w[2]) };
    std::sort(sorted, sorted + 3);
    EXPECT_NEAR(1.0, sorted[0], 10 * kEps);
    EXPECT_NEAR(3.0, sorted[1], 10 * kEps);
    EXPECT_NEAR(5.0, sorted[2], 10 * kEps);
    EXPECT_EQ(99, a[3]);
}

/* Degenerate inputs: diagonal, repeated eigenvalues, zero, huge, tiny */
TEST(LinAlg, EigenSym3Degenerate) {
    const trical_wide_t cases[5][9] = {
        { 3, 0, 0, 0, 1, 0, 0, 0, 2 },
        { 2, 0, 0, 0, 2, 0, 0, 0, 2 },
        { 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 4, 2, 0, 2, 4, 0, 0, 0, 2 }
    };
    const double scales[3] = { 1.0, 1e15, 1e-15 };

    for (int c = 0; c < 5; c++) {
        for (int s = 0; s < 3; s++) {
            trical_wide_t a[9], w[3], v[9];
            double trace = 0.0, sum = 0.0;

            for (int i = 0; i < 9; i++) {
                a[i] = static_cast<trical_wide_t>(cases[c][i] * scales[s]);
            }
            trical_eigen_sym3(a, w, v);

            for (int i = 0; i < 3; i++) {
                EXPECT_TRUE(std::isfinite(static_cast<double>(w[i])));
                trace += static_cast<double>(a[i * 4]);
                sum += static_cast<double>(w[i]);
            }
            for (int i = 0; i < 9; i++) {
                EXPECT_TRUE(std::isfinite(static_cast<double>(v[i])));
            }
            EXPECT_NEAR(trace, sum, 50 * kEps * scales[s] + 1e-300);
        }
    }
}

TEST(LinAlg, Det3) {
    const trical_wide_t a[9] = { 2, 0, 1, 1, 3, 2, 1, 1, 1 };
    const trical_wide_t identity[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

    /* 2 (3 - 2) - 0 + 1 (1 - 3) = 0 */
    EXPECT_NEAR(0.0, trical_det3(a), kEps);
    EXPECT_NEAR(1.0, trical_det3(identity), kEps);
}
