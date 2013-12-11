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
#include <cstdint>
#include <cstdlib>
#include <cstdbool>
#include <cassert>
#include <cmath>
#include <cfloat>

/* C++ complains about the C99 'restrict' qualifier. Just ignore it. */
#define restrict

#include "TRICAL.h"
#include "3dmath.h"

/*
Test LLT decomposition of a square positive definite matrix via the Cholesky
method.

The expected output is a lower-triangle matrix in column-major format.
*/
TEST(MatrixMath, CholeskyLLT) {
    float m1[16] = {
        18, 22,  54,  42,
        22, 70,  86,  62,
        54, 86, 174, 134,
        42, 62, 134, 106
    };
    float expected[16] = {
        4.24264, 5.18545, 12.72792, 9.8994951,
        0, 6.5659051, 3.0460384, 1.6245539,
        0, 0, 1.6497422, 1.8497111,
        0, 0, 0, 1.3926213
    };

    matrix_cholesky_decomp_scale_f(4, m1, m1, 1.0);
    EXPECT_NEAR(expected[0], m1[0], 1e-5);
    EXPECT_NEAR(expected[1], m1[1], 1e-5);
    EXPECT_NEAR(expected[2], m1[2], 1e-5);
    EXPECT_NEAR(expected[3], m1[3], 1e-5);
    EXPECT_NEAR(expected[5], m1[5], 1e-5);
    EXPECT_NEAR(expected[6], m1[6], 1e-5);
    EXPECT_NEAR(expected[7], m1[7], 1e-5);
    EXPECT_NEAR(expected[10], m1[10], 1e-5);
    EXPECT_NEAR(expected[11], m1[11], 1e-5);
    EXPECT_NEAR(expected[15], m1[15], 1e-5);
}

/*
Test multiplication of square and rectangular matricies
*/
TEST(MatrixMath, Multiply) {
    float a[16] = {
        1, 1, 1, 1,
        2, 4, 8, 16,
        3, 9, 27, 81,
        4, 16, 64, 256
    };
    float b[16] = {
        4.0, -3.0, 4.0/3.0, -1.0/4.0,
        -13.0/3.0, 19.0/4.0, -7.0/3.0, 11.0/24.0,
        3.0/2.0, -2.0, 7.0/6.0, -1.0/4.0,
        -1.0/6.0, 1.0/4.0, -1.0/6.0, 1.0/24.0
    };
    float c[16];

    matrix_multiply_f(c, a, b, 4, 4, 4, 4, 1.0);
    EXPECT_NEAR(1.0, c[0], 1e-5);
    EXPECT_NEAR(0.0, c[1], 1e-5);
    EXPECT_NEAR(0.0, c[2], 1e-5);
    EXPECT_NEAR(0.0, c[3], 1e-5);
    EXPECT_NEAR(0.0, c[4], 1e-5);
    EXPECT_NEAR(1.0, c[5], 1e-5);
    EXPECT_NEAR(0.0, c[6], 1e-5);
    EXPECT_NEAR(0.0, c[7], 1e-5);
    EXPECT_NEAR(0.0, c[8], 1e-5);
    EXPECT_NEAR(0.0, c[9], 1e-5);
    EXPECT_NEAR(1.0, c[10], 1e-5);
    EXPECT_NEAR(0.0, c[11], 1e-5);
    EXPECT_NEAR(0.0, c[12], 1e-5);
    EXPECT_NEAR(0.0, c[13], 1e-5);
    EXPECT_NEAR(0.0, c[14], 1e-5);
    EXPECT_NEAR(1.0, c[15], 1e-5);

    float d[4] = { 1, 2, 3, 4}, e[4] = { 5, 6, 7, 8 };

    matrix_multiply_f(c, d, e, 1, 4, 4, 1, 1.0);
    EXPECT_NEAR(5.0, c[0], 1e-5);
    EXPECT_NEAR(10.0, c[1], 1e-5);
    EXPECT_NEAR(15.0, c[2], 1e-5);
    EXPECT_NEAR(20.0, c[3], 1e-5);
    EXPECT_NEAR(6.0, c[4], 1e-5);
    EXPECT_NEAR(12.0, c[5], 1e-5);
    EXPECT_NEAR(18.0, c[6], 1e-5);
    EXPECT_NEAR(24.0, c[7], 1e-5);
    EXPECT_NEAR(7.0, c[8], 1e-5);
    EXPECT_NEAR(14.0, c[9], 1e-5);
    EXPECT_NEAR(21.0, c[10], 1e-5);
    EXPECT_NEAR(28.0, c[11], 1e-5);
    EXPECT_NEAR(8.0, c[12], 1e-5);
    EXPECT_NEAR(16.0, c[13], 1e-5);
    EXPECT_NEAR(24.0, c[14], 1e-5);
    EXPECT_NEAR(32.0, c[15], 1e-5);
}
