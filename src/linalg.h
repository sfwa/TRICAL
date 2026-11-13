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

#ifndef TRICAL_LINALG_H
#define TRICAL_LINALG_H

/*
Small dense linear algebra kernels in the "wide" floating-point type, for
the handful of fixed-size problems TRICAL solves. All matrices are row-major.
*/

#include <math.h>

#include "TRICAL.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef TRICAL_WIDE_TYPE trical_wide_t;

/*
Square root and absolute value in the wide type. sizeof is evaluated at
compile time, so only one branch of each survives.
*/
static inline trical_wide_t trical_sqrt(trical_wide_t x) {
    return (sizeof(trical_wide_t) == sizeof(float)) ?
        (trical_wide_t)sqrtf((float)x) : (trical_wide_t)sqrt((double)x);
}

static inline trical_wide_t trical_abs(trical_wide_t x) {
    return x < 0 ? -x : x;
}

/*
trical_cholesky
In-place Cholesky factorisation of the symmetric positive definite `n` x `n`
matrix `a`. Only the lower triangle of `a` is read; on success it holds L
such that L L^T = a (the strict upper triangle is left as it was).

Returns 0 on success, and -1 if `a` is not (numerically) positive definite,
in which case its contents are unspecified.
*/
static inline int trical_cholesky(trical_wide_t *a, unsigned int n) {
    unsigned int i, j, k;

    for (j = 0; j < n; j++) {
        trical_wide_t d = a[j * n + j], r;

        for (k = 0; k < j; k++) {
            d -= a[j * n + k] * a[j * n + k];
        }
        /* written so that a NaN pivot is rejected as well */
        if (!(d > 0)) {
            return -1;
        }
        d = trical_sqrt(d);
        a[j * n + j] = d;
        r = (trical_wide_t)1 / d;

        for (i = j + 1; i < n; i++) {
            trical_wide_t s = a[i * n + j];

            for (k = 0; k < j; k++) {
                s -= a[i * n + k] * a[j * n + k];
            }
            a[i * n + j] = s * r;
        }
    }

    return 0;
}

/*
trical_cholesky_solve
Solves (L L^T) x = b in place, given the factor L produced by trical_cholesky
in the lower triangle of `l`. On entry `x` holds b; on return it holds x.
*/
static inline void trical_cholesky_solve(const trical_wide_t *l,
unsigned int n, trical_wide_t *x) {
    unsigned int i, k;

    /* forward substitution: L y = b */
    for (i = 0; i < n; i++) {
        trical_wide_t s = x[i];

        for (k = 0; k < i; k++) {
            s -= l[i * n + k] * x[k];
        }
        x[i] = s / l[i * n + i];
    }

    /* back substitution: L^T x = y */
    for (i = n; i-- > 0;) {
        trical_wide_t s = x[i];

        for (k = i + 1; k < n; k++) {
            s -= l[k * n + i] * x[k];
        }
        x[i] = s / l[i * n + i];
    }
}

/*
trical_eigen_sym3
Eigen-decomposition of the symmetric 3x3 matrix `a` by cyclic Jacobi
rotations: on return `w` holds the eigenvalues (in no particular order) and
the columns of `v` the corresponding orthonormal eigenvectors, so that
a = v diag(w) v^T. Only the upper triangle of `a` is read; `a` is not
modified.

The iteration count is bounded, and the result is finite whenever the input
is.
*/
static inline void trical_eigen_sym3(const trical_wide_t a[9],
trical_wide_t w[3], trical_wide_t v[9]) {
    /* the three off-diagonal pairs, and the remaining index for each */
    static const unsigned char pairs[3][3] = {
        { 0, 1, 2 }, { 0, 2, 1 }, { 1, 2, 0 }
    };
    trical_wide_t m[9];
    unsigned int sweep, i;

    m[0] = a[0]; m[1] = a[1]; m[2] = a[2];
    m[3] = a[1]; m[4] = a[4]; m[5] = a[5];
    m[6] = a[2]; m[7] = a[5]; m[8] = a[8];

    v[0] = 1; v[1] = 0; v[2] = 0;
    v[3] = 0; v[4] = 1; v[5] = 0;
    v[6] = 0; v[7] = 0; v[8] = 1;

    /*
    Cyclic Jacobi converges quadratically; for a 3x3 matrix 5 or 6 sweeps
    reach working precision, in single or double precision alike.
    */
    for (sweep = 0; sweep < 12; sweep++) {
        trical_wide_t off = trical_abs(m[1]) + trical_abs(m[2]) +
                            trical_abs(m[5]);
        trical_wide_t diag = trical_abs(m[0]) + trical_abs(m[4]) +
                             trical_abs(m[8]);

        /* done once the off-diagonal part no longer registers */
        if (diag + off == diag) {
            break;
        }

        for (i = 0; i < 3; i++) {
            unsigned int p = pairs[i][0], q = pairs[i][1], r = pairs[i][2];
            trical_wide_t apq = m[p * 3 + q], theta, t, c, s, tmp;
            unsigned int k;

            if (apq == 0) {
                continue;
            }

            /*
            Rotation angle that zeroes m[p][q]: t = tan(angle) is the
            smaller root of t^2 + 2 t theta - 1 = 0.
            */
            theta = (m[q * 3 + q] - m[p * 3 + p]) / (apq + apq);
            t = (trical_wide_t)1 /
                (trical_abs(theta) + trical_sqrt(theta * theta + 1));
            if (theta < 0) {
                t = -t;
            }
            c = (trical_wide_t)1 / trical_sqrt(t * t + 1);
            s = t * c;

            m[p * 3 + p] -= t * apq;
            m[q * 3 + q] += t * apq;
            m[p * 3 + q] = 0;
            m[q * 3 + p] = 0;

            tmp = m[r * 3 + p];
            m[r * 3 + p] = c * tmp - s * m[r * 3 + q];
            m[r * 3 + q] = s * tmp + c * m[r * 3 + q];
            m[p * 3 + r] = m[r * 3 + p];
            m[q * 3 + r] = m[r * 3 + q];

            for (k = 0; k < 3; k++) {
                tmp = v[k * 3 + p];
                v[k * 3 + p] = c * tmp - s * v[k * 3 + q];
                v[k * 3 + q] = s * tmp + c * v[k * 3 + q];
            }
        }
    }

    w[0] = m[0];
    w[1] = m[4];
    w[2] = m[8];
}

/* 3x3 determinant */
static inline trical_wide_t trical_det3(const trical_wide_t a[9]) {
    return a[0] * (a[4] * a[8] - a[5] * a[7]) -
           a[1] * (a[3] * a[8] - a[5] * a[6]) +
           a[2] * (a[3] * a[7] - a[4] * a[6]);
}

#ifdef __cplusplus
}
#endif

#endif
