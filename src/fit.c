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
#include <math.h>
#include <stdint.h>
#include <string.h>

#include "TRICAL.h"
#include "fit.h"
#include "linalg.h"

/*
A bit about the formulation in this file. The main references are
[1]: Kanatani, "Overviews of Optimization Techniques for Geometric
     Estimation", Mem. Fac. Eng. Okayama Univ. 47 (2013)
[2]: Markovsky, Kukush & Van Huffel, "Consistent least squares fitting of
     ellipsoids", Numer. Math. 98 (2004)
[3]: Taubin, "Estimation of planar curves, surfaces, and nonplanar space
     curves defined by implicit equations...", IEEE PAMI 13 (1991)
[4]: Vasconcelos et al., "Geometric approach to strapdown magnetometer
     calibration in sensor frame", IEEE TAES 47 (2011)
[5]: Markley, "Attitude determination using vector observations and the
     singular value decomposition", J. Astronaut. Sci. 36 (1988)

A sensor reading B is related to the field H in the sensor frame by
B = (I + D)^-1 H + b + noise. |H| is known (the field norm), so the readings
lie on an ellipsoid, which determines b and the symmetric positive definite
part of (I + D) -- everything except a rotation ([4]).

Working in normalised coordinates x = (B - origin) / scale, the ellipsoid is
the zero set of a quadric

    F(x) = x' G x + 2 q' x + q0 = xi(x) . theta
    xi(x) = (xx, yy, zz, 2yz, 2xz, 2xy, 2x, 2y, 2z, 1)
    theta = (G11, G22, G33, G23, G13, G12, q1, q2, q3, q0)

which is linear in theta. Everything a least-squares fit of theta needs is
the moment matrix M = sum xi xi', whose entries are sums of monomials of x of
degree <= 4: 35 running sums, however many readings there are. That is all
TRICAL stores, so adding a reading costs a few dozen multiply-adds, and the
sums do not depend on the order of the readings or on any estimate that was
current when they arrived.

Minimising theta' M theta alone is biased, because the noise in x is squared
by xi: E[M] = M_true + s^2 N + O(s^4), where s is the noise standard
deviation and N is known in terms of the data ([1], sec. 4.2-4.4). So
instead of the smallest eigenvector of M we want the one of the generalised
problem

    M theta = lambda N theta

which removes the bias (correcting least squares in this way is what makes
it consistent, [2]); the eigenvalue lambda is then an estimate of s^2. With
N = sum J J' (J the Jacobian of xi) this is Taubin's method [3], which
normalises the algebraic distance F by its gradient to approximate the
geometric distance: theta' M theta / theta' N theta is the mean squared
distance of the readings from the ellipsoid. Adding the term 2 S[xi e'],
which accounts for the bias of xi itself (E[xi] = xi_true + s^2 e), gives the
leading part of the "hyper-accurate" normalisation of [1], correct to second
order in the noise. Its accuracy matches a maximum-likelihood fit (see
test/test_accuracy.cpp).

The problem is solved in the form in which the "renormalisation" method was
first given ([1], sec. 4.2 and its footnote): theta is the eigenvector for
the smallest eigenvalue of M - lambda N, where lambda is the Rayleigh
quotient theta' M theta / theta' N theta of the previous solution. That
iteration converges quickly, and since the fit is refreshed as readings
arrive, each refit normally takes just one step of it.

Where the readings do not determine the ellipsoid (few readings, a sensor at
rest, motion about one axis only) the problem is ill-posed, so a prior holds
theta to the initial estimate: a positive semi-definite matrix w Q is added
to M - lambda N, which penalises differences from the initial scale factors
and, less firmly, from the initial bias (see _prior_t below). That is an
ordinary eigenproblem still. Its weight w is in proportion to s^2, for which
the estimate above is used, so the prior gives way to the readings at the
rate that their noise warrants -- and it is much the greater while the
readings do not contradict the initial estimate, since noise would otherwise
push about whatever they leave undetermined (see _regularisation).

The rotation that the ellipsoid leaves open is taken to be the identity (so
that D is symmetric) unless reference directions are supplied, in which case
it is the solution of Wahba's problem [5] for the pairs (reference direction,
calibrated reading). That only needs the sums of f x' and of f, where f is
the reference direction, so it is likewise independent of the order of the
readings and always consistent with the current ellipsoid.
*/

#define TRICAL_WIDE_EPSILON ((sizeof(trical_wide_t) == sizeof(float)) ? \
                             (trical_wide_t)FLT_EPSILON : \
                             (trical_wide_t)DBL_EPSILON)

/* Tables generated by tools/gen_tables.py */

/* exponents (a, b, c) of monomial k = x^a y^b z^c */
static const unsigned char TRICAL_MONO_EXP[TRICAL_NUM_MOMENTS][3] = {
    {  0,  0,  0 },
    {  1,  0,  0 },
    {  0,  1,  0 },
    {  0,  0,  1 },
    {  2,  0,  0 },
    {  1,  1,  0 },
    {  1,  0,  1 },
    {  0,  2,  0 },
    {  0,  1,  1 },
    {  0,  0,  2 },
    {  3,  0,  0 },
    {  2,  1,  0 },
    {  2,  0,  1 },
    {  1,  2,  0 },
    {  1,  1,  1 },
    {  1,  0,  2 },
    {  0,  3,  0 },
    {  0,  2,  1 },
    {  0,  1,  2 },
    {  0,  0,  3 },
    {  4,  0,  0 },
    {  3,  1,  0 },
    {  3,  0,  1 },
    {  2,  2,  0 },
    {  2,  1,  1 },
    {  2,  0,  2 },
    {  1,  3,  0 },
    {  1,  2,  1 },
    {  1,  1,  2 },
    {  1,  0,  3 },
    {  0,  4,  0 },
    {  0,  3,  1 },
    {  0,  2,  2 },
    {  0,  1,  3 },
    {  0,  0,  4 },
};

/* index of the monomial that is the product of quadric basis terms i, j */
static const unsigned char TRICAL_M_INDEX[10][10] = {
    { 20, 23, 25, 24, 22, 21, 10, 11, 12,  4 },
    { 23, 30, 32, 31, 27, 26, 13, 16, 17,  7 },
    { 25, 32, 34, 33, 29, 28, 15, 18, 19,  9 },
    { 24, 31, 33, 32, 28, 27, 14, 17, 18,  8 },
    { 22, 27, 29, 28, 25, 24, 12, 14, 15,  6 },
    { 21, 26, 28, 27, 24, 23, 11, 13, 14,  5 },
    { 10, 13, 15, 14, 12, 11,  4,  5,  6,  1 },
    { 11, 16, 18, 17, 14, 13,  5,  7,  8,  2 },
    { 12, 17, 19, 18, 15, 14,  6,  8,  9,  3 },
    {  4,  7,  9,  8,  6,  5,  1,  2,  3,  0 },
};

/* integer coefficient of quadric basis term i */
static const unsigned char TRICAL_XI_COEF[10] = {
    1, 1, 1, 2, 2, 2, 2, 2, 2, 1
};

/* index of monomial x^a y^b z^c for a + b + c <= 4 (255 otherwise) */
static const unsigned char TRICAL_MONO_INDEX[5][5][5] = {
    {
        {   0,   3,   9,  19,  34 },
        {   2,   8,  18,  33, 255 },
        {   7,  17,  32, 255, 255 },
        {  16,  31, 255, 255, 255 },
        {  30, 255, 255, 255, 255 },
    },
    {
        {   1,   6,  15,  29, 255 },
        {   5,  14,  28, 255, 255 },
        {  13,  27, 255, 255, 255 },
        {  26, 255, 255, 255, 255 },
        { 255, 255, 255, 255, 255 },
    },
    {
        {   4,  12,  25, 255, 255 },
        {  11,  24, 255, 255, 255 },
        {  23, 255, 255, 255, 255 },
        { 255, 255, 255, 255, 255 },
        { 255, 255, 255, 255, 255 },
    },
    {
        {  10,  22, 255, 255, 255 },
        {  21, 255, 255, 255, 255 },
        { 255, 255, 255, 255, 255 },
        { 255, 255, 255, 255, 255 },
        { 255, 255, 255, 255, 255 },
    },
    {
        {  20, 255, 255, 255, 255 },
        { 255, 255, 255, 255, 255 },
        { 255, 255, 255, 255, 255 },
        { 255, 255, 255, 255, 255 },
        { 255, 255, 255, 255, 255 },
    },
};

/*
True if `x` is neither infinite nor NaN. These work on the bit pattern, so
that they cannot be optimised away by compilers told to assume finite
arithmetic.
*/
static int _float_isfinite(float x) {
    uint32_t bits;

    memcpy(&bits, &x, sizeof(bits));
    return (bits & 0x7f800000u) != 0x7f800000u;
}

static int _wide_isfinite(trical_wide_t x) {
    if (sizeof(trical_wide_t) == sizeof(float)) {
        return _float_isfinite((float)x);
    } else {
        const double value = (double)x;
        uint64_t bits;

        memcpy(&bits, &value, sizeof(bits));
        return (bits >> 52 & 0x7ffu) != 0x7ffu;
    }
}

/* True if `x` is a number that can be converted to single precision. */
static int _wide_isfloat(trical_wide_t x) {
    return _wide_isfinite(x) && x <= (trical_wide_t)FLT_MAX &&
           x >= -(trical_wide_t)FLT_MAX;
}

static trical_wide_t _dot10(const trical_wide_t a[10],
const trical_wide_t b[10]) {
    trical_wide_t s = 0;
    unsigned int i;

    for (i = 0; i < 10; i++) {
        s += a[i] * b[i];
    }

    return s;
}

/* c = a b for 3x3 matrices; c must not alias a or b */
static void _mul33(const trical_wide_t a[9], const trical_wide_t b[9],
trical_wide_t c[9]) {
    unsigned int i, j;

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            c[i * 3 + j] = a[i * 3 + 0] * b[0 * 3 + j] +
                           a[i * 3 + 1] * b[1 * 3 + j] +
                           a[i * 3 + 2] * b[2 * 3 + j];
        }
    }
}

/* out = v diag(d) v^T for the 3x3 matrix v; exactly symmetric */
static void _similarity3(const trical_wide_t v[9], const trical_wide_t d[3],
trical_wide_t out[9]) {
    unsigned int i, j;

    for (i = 0; i < 3; i++) {
        for (j = i; j < 3; j++) {
            out[i * 3 + j] = v[i * 3 + 0] * d[0] * v[j * 3 + 0] +
                             v[i * 3 + 1] * d[1] * v[j * 3 + 1] +
                             v[i * 3 + 2] * d[2] * v[j * 3 + 2];
            out[j * 3 + i] = out[i * 3 + j];
        }
    }
}

/*
_trical_moments_accumulate
Adds the 35 monomials of degree <= 4 of the normalised measurement `x` to
`sums`. The products are formed in single precision.
*/
void _trical_moments_accumulate(trical_wide_t sums[TRICAL_NUM_MOMENTS],
const float v[3]) {
    const float x = v[0], y = v[1], z = v[2];
    const float xx = x * x, xy = x * y, xz = x * z;
    const float yy = y * y, yz = y * z, zz = z * z;
    trical_wide_t *restrict s = sums;

    s[0] += (trical_wide_t)1;

    s[1] += (trical_wide_t)x;
    s[2] += (trical_wide_t)y;
    s[3] += (trical_wide_t)z;

    s[4] += (trical_wide_t)xx;
    s[5] += (trical_wide_t)xy;
    s[6] += (trical_wide_t)xz;
    s[7] += (trical_wide_t)yy;
    s[8] += (trical_wide_t)yz;
    s[9] += (trical_wide_t)zz;

    s[10] += (trical_wide_t)(xx * x);
    s[11] += (trical_wide_t)(xx * y);
    s[12] += (trical_wide_t)(xx * z);
    s[13] += (trical_wide_t)(x * yy);
    s[14] += (trical_wide_t)(xy * z);
    s[15] += (trical_wide_t)(x * zz);
    s[16] += (trical_wide_t)(yy * y);
    s[17] += (trical_wide_t)(yy * z);
    s[18] += (trical_wide_t)(y * zz);
    s[19] += (trical_wide_t)(zz * z);

    s[20] += (trical_wide_t)(xx * xx);
    s[21] += (trical_wide_t)(xx * xy);
    s[22] += (trical_wide_t)(xx * xz);
    s[23] += (trical_wide_t)(xx * yy);
    s[24] += (trical_wide_t)(xx * yz);
    s[25] += (trical_wide_t)(xx * zz);
    s[26] += (trical_wide_t)(xy * yy);
    s[27] += (trical_wide_t)(xy * yz);
    s[28] += (trical_wide_t)(xy * zz);
    s[29] += (trical_wide_t)(xz * zz);
    s[30] += (trical_wide_t)(yy * yy);
    s[31] += (trical_wide_t)(yy * yz);
    s[32] += (trical_wide_t)(yy * zz);
    s[33] += (trical_wide_t)(yz * zz);
    s[34] += (trical_wide_t)(zz * zz);
}

/*
_trical_moments_shift
Converts sums of monomials of x to sums of monomials of (x - delta), by
expanding each (x - delta)^a (y - delta)^b (z - delta)^c binomially.
*/
void _trical_moments_shift(trical_wide_t sums[TRICAL_NUM_MOMENTS],
const trical_wide_t delta[3]) {
    static const unsigned char binomial[5][5] = {
        { 1, 0, 0, 0, 0 },
        { 1, 1, 0, 0, 0 },
        { 1, 2, 1, 0, 0 },
        { 1, 3, 3, 1, 0 },
        { 1, 4, 6, 4, 1 }
    };
    trical_wide_t power[3][5], out[TRICAL_NUM_MOMENTS];
    unsigned int i, k;

    assert(sums && delta);

    /* power[i][k] = (-delta[i])^k */
    for (i = 0; i < 3; i++) {
        power[i][0] = 1;
        for (k = 1; k < 5; k++) {
            power[i][k] = -power[i][k - 1] * delta[i];
        }
    }

    for (k = 0; k < TRICAL_NUM_MOMENTS; k++) {
        const unsigned int a = TRICAL_MONO_EXP[k][0];
        const unsigned int b = TRICAL_MONO_EXP[k][1];
        const unsigned int c = TRICAL_MONO_EXP[k][2];
        trical_wide_t acc = 0;
        unsigned int p, q, r;

        for (p = 0; p <= a; p++) {
            for (q = 0; q <= b; q++) {
                for (r = 0; r <= c; r++) {
                    const unsigned int count = (unsigned int)binomial[a][p] *
                        binomial[b][q] * binomial[c][r];

                    acc += (trical_wide_t)count * power[0][a - p] *
                           power[1][b - q] * power[2][c - r] *
                           sums[TRICAL_MONO_INDEX[p][q][r]];
                }
            }
        }
        out[k] = acc;
    }

    memcpy(sums, out, sizeof(out));
}

/*
_trical_basis_mean
Writes to `xi` the mean of the quadric basis
xi = (xx, yy, zz, 2yz, 2xz, 2xy, 2x, 2y, 2z, 1) over the measurements, from
the monomial sums.
*/
void _trical_basis_mean(const trical_wide_t sums[TRICAL_NUM_MOMENTS],
trical_wide_t xi[10]) {
    const trical_wide_t scale = (trical_wide_t)1 / sums[0];
    unsigned int i;

    assert(sums && xi);

    for (i = 0; i < 10; i++) {
        xi[i] = (trical_wide_t)TRICAL_XI_COEF[i] *
                sums[TRICAL_M_INDEX[i][9]] * scale;
    }
}

/*
_trical_moment_product
y = M theta, where M = E[xi xi^T] is the mean moment matrix of the quadric
basis. Its entries are monomial sums (with small integer coefficients), so
the matrix itself is never formed.
*/
void _trical_moment_product(const trical_wide_t sums[TRICAL_NUM_MOMENTS],
const trical_wide_t theta[10], trical_wide_t y[10]) {
    const trical_wide_t scale = (trical_wide_t)1 / sums[0];
    trical_wide_t scaled[10];
    unsigned int i, j;

    assert(sums && theta && y);

    for (j = 0; j < 10; j++) {
        scaled[j] = (trical_wide_t)TRICAL_XI_COEF[j] * theta[j];
    }
    for (i = 0; i < 10; i++) {
        const unsigned char *index = TRICAL_M_INDEX[i];
        trical_wide_t s = 0;

        for (j = 0; j < 10; j++) {
            s += sums[index[j]] * scaled[j];
        }
        y[i] = (trical_wide_t)TRICAL_XI_COEF[i] * s * scale;
    }
}

/*
The Taubin normalisation matrix N_T = E[J J^T], where J is the Jacobian of
xi with respect to x. The rows of J are
    xx: (2x, 0, 0)    yy: (0, 2y, 0)    zz: (0, 0, 2z)
    2yz: (0, 2z, 2y)  2xz: (2z, 0, 2x)  2xy: (2y, 2x, 0)
    2x: (2, 0, 0)     2y: (0, 2, 0)     2z: (0, 0, 2)     1: (0, 0, 0)
so the entries of J J^T are (multiples of) first and second moments of x,
i.e. of elements of the mean of xi. Each row here is one term of the upper
triangle: { row, column, element of the mean of xi, coefficient }.
*/
static const unsigned char TRICAL_N_TERMS[30][4] = {
    { 0, 0, 0, 4 }, { 0, 4, 4, 2 }, { 0, 5, 5, 2 }, { 0, 6, 6, 2 },
    { 1, 1, 1, 4 }, { 1, 3, 3, 2 }, { 1, 5, 5, 2 }, { 1, 7, 7, 2 },
    { 2, 2, 2, 4 }, { 2, 3, 3, 2 }, { 2, 4, 4, 2 }, { 2, 8, 8, 2 },
    { 3, 3, 1, 4 }, { 3, 3, 2, 4 }, { 3, 4, 5, 2 }, { 3, 5, 4, 2 },
    { 3, 7, 8, 2 }, { 3, 8, 7, 2 },
    { 4, 4, 0, 4 }, { 4, 4, 2, 4 }, { 4, 5, 3, 2 }, { 4, 6, 8, 2 },
    { 4, 8, 6, 2 },
    { 5, 5, 0, 4 }, { 5, 5, 1, 4 }, { 5, 6, 7, 2 }, { 5, 7, 6, 2 },
    { 6, 6, 9, 4 }, { 7, 7, 9, 4 }, { 8, 8, 9, 4 }
};

/*
_trical_normalisation_product
y = N theta, where N is the Taubin normalisation matrix N_T for the
measurements whose mean quadric basis is `xi` -- plus, if `hyper` is
non-zero, the second-order ("hyper") correction xi e^T + e xi^T, with
e = (1, 1, 1, 0, ..., 0).
*/
void _trical_normalisation_product(const trical_wide_t xi[10],
const trical_wide_t theta[10], int hyper, trical_wide_t y[10]) {
    unsigned int i, k;

    assert(xi && theta && y);

    for (i = 0; i < 10; i++) {
        y[i] = 0;
    }
    for (k = 0; k < 30; k++) {
        const unsigned int row = TRICAL_N_TERMS[k][0],
                           column = TRICAL_N_TERMS[k][1];
        const trical_wide_t value = (trical_wide_t)TRICAL_N_TERMS[k][3] *
                                    xi[TRICAL_N_TERMS[k][2]];

        y[row] += value * theta[column];
        if (row != column) {
            y[column] += value * theta[row];
        }
    }

    if (hyper) {
        const trical_wide_t e_theta = theta[0] + theta[1] + theta[2];
        const trical_wide_t xi_theta = _dot10(xi, theta);

        for (i = 0; i < 10; i++) {
            y[i] += xi[i] * e_theta;
        }
        for (i = 0; i < 3; i++) {
            y[i] += xi_theta;
        }
    }
}

/* Scales `quadric` to unit Euclidean norm. */
static void _quadric_normalise(trical_wide_t quadric[10]) {
    const trical_wide_t norm = trical_sqrt(_dot10(quadric, quadric));
    unsigned int i;

    if (norm > 0) {
        for (i = 0; i < 10; i++) {
            quadric[i] /= norm;
        }
    }
}

/*
_trical_quadric_from_calibration
Writes to `quadric` the unit-norm coefficients of the ellipsoid
|A (x - c)| = 1, where `ata` is A^T A:
    F(x) = (x - c)' G (x - c) - 1, i.e. q = -G c, q0 = c' G c - 1
*/
void _trical_quadric_from_calibration(const trical_wide_t ata[9],
const trical_wide_t centre[3], trical_wide_t quadric[10]) {
    trical_wide_t gc[3];
    unsigned int i;

    assert(ata && centre && quadric);

    for (i = 0; i < 3; i++) {
        gc[i] = ata[i * 3 + 0] * centre[0] + ata[i * 3 + 1] * centre[1] +
                ata[i * 3 + 2] * centre[2];
    }

    quadric[0] = ata[0 * 3 + 0];
    quadric[1] = ata[1 * 3 + 1];
    quadric[2] = ata[2 * 3 + 2];
    quadric[3] = ata[1 * 3 + 2];
    quadric[4] = ata[0 * 3 + 2];
    quadric[5] = ata[0 * 3 + 1];
    quadric[6] = -gc[0];
    quadric[7] = -gc[1];
    quadric[8] = -gc[2];
    quadric[9] = centre[0] * gc[0] + centre[1] * gc[1] + centre[2] * gc[2] -
                 1;

    _quadric_normalise(quadric);
}

/*
_trical_quadric_to_calibration
Finds the centre c and symmetric positive definite A such that the zero set
of `quadric` is |A (x - c)| = 1:
    c = -G^-1 q,  A = (G / rho)^(1/2)  with rho = q' G^-1 q - q0

Returns 0 on success, or -1 if `quadric` is not an ellipsoid of acceptable
shape and size. `a_inverse` may be NULL.
*/
int _trical_quadric_to_calibration(const trical_wide_t quadric[10],
trical_wide_t centre[3], trical_wide_t a[9], trical_wide_t a_inverse[9]) {
    const trical_wide_t ratio2 = (trical_wide_t)(TRICAL_AXIS_RATIO_MAX *
                                                 TRICAL_AXIS_RATIO_MAX);
    trical_wide_t g[9], q[3], q0, w[3], v[9], d[3], vtq[3], rho, wmin, wmax;
    trical_wide_t sign = 1;
    unsigned int i;

    assert(quadric && centre && a);

    /* a quadric is only defined up to scale; make the trace of G positive */
    if (quadric[0] + quadric[1] + quadric[2] < 0) {
        sign = -1;
    }
    g[0] = sign * quadric[0];
    g[4] = sign * quadric[1];
    g[8] = sign * quadric[2];
    g[5] = g[7] = sign * quadric[3];
    g[2] = g[6] = sign * quadric[4];
    g[1] = g[3] = sign * quadric[5];
    q[0] = sign * quadric[6];
    q[1] = sign * quadric[7];
    q[2] = sign * quadric[8];
    q0 = sign * quadric[9];

    trical_eigen_sym3(g, w, v);

    wmin = w[0] < w[1] ? (w[0] < w[2] ? w[0] : w[2]) :
                         (w[1] < w[2] ? w[1] : w[2]);
    wmax = w[0] > w[1] ? (w[0] > w[2] ? w[0] : w[2]) :
                         (w[1] > w[2] ? w[1] : w[2]);

    /* an ellipsoid needs G positive definite; reject extreme shapes too */
    if (!(wmin > 0) || !(wmax <= ratio2 * wmin)) {
        return -1;
    }

    /* c = -G^-1 q and rho = q' G^-1 q - q0, via the eigenvectors of G */
    rho = -q0;
    for (i = 0; i < 3; i++) {
        vtq[i] = v[0 * 3 + i] * q[0] + v[1 * 3 + i] * q[1] +
                 v[2 * 3 + i] * q[2];
        rho += vtq[i] * vtq[i] / w[i];
    }
    if (!(rho > 0) || !_wide_isfinite(rho)) {
        return -1;
    }

    for (i = 0; i < 3; i++) {
        centre[i] = -(v[i * 3 + 0] * vtq[0] / w[0] +
                      v[i * 3 + 1] * vtq[1] / w[1] +
                      v[i * 3 + 2] * vtq[2] / w[2]);
        if (!_wide_isfinite(centre[i])) {
            return -1;
        }
    }

    /*
    The principal scale factors, which are the reciprocals of the semi-axes.
    The normalised coordinates are scaled so that these are about 1 for the
    initial estimate; a fit that is wildly different (e.g. a nearly flat
    "ellipsoid" through a small patch of measurements, or a tiny one around
    a few of them) is not credible.
    */
    for (i = 0; i < 3; i++) {
        d[i] = trical_sqrt(w[i] / rho);
        if (!_wide_isfinite(d[i]) ||
                !(d[i] * (trical_wide_t)TRICAL_SIZE_MAX > 1) ||
                !(d[i] * (trical_wide_t)TRICAL_SIZE_MIN < 1)) {
            return -1;
        }
    }
    _similarity3(v, d, a);

    if (a_inverse) {
        for (i = 0; i < 3; i++) {
            d[i] = (trical_wide_t)1 / d[i];
        }
        _similarity3(v, d, a_inverse);
    }

    return 0;
}

/*
_trical_quadric_translate
Re-expresses `quadric` in coordinates x' = x - delta, and normalises it:
    F = (x' + delta)' G (x' + delta) + 2 q' (x' + delta) + q0
so G is unchanged, q' = q + G delta, q0' = q0 + 2 q' delta + delta' G delta.
*/
void _trical_quadric_translate(trical_wide_t quadric[10],
const trical_wide_t delta[3]) {
    trical_wide_t gd[3];

    assert(quadric && delta);

    gd[0] = quadric[0] * delta[0] + quadric[5] * delta[1] +
            quadric[4] * delta[2];
    gd[1] = quadric[5] * delta[0] + quadric[1] * delta[1] +
            quadric[3] * delta[2];
    gd[2] = quadric[4] * delta[0] + quadric[3] * delta[1] +
            quadric[2] * delta[2];

    quadric[9] += 2 * (quadric[6] * delta[0] + quadric[7] * delta[1] +
                       quadric[8] * delta[2]) +
                  delta[0] * gd[0] + delta[1] * gd[1] + delta[2] * gd[2];
    quadric[6] += gd[0];
    quadric[7] += gd[1];
    quadric[8] += gd[2];

    _quadric_normalise(quadric);
}

/*
_trical_polar_rotation
Writes to `q` the rotation matrix closest to `k`: the solution of Wahba's
problem for the attitude profile matrix `k`.

With the singular value decomposition k = U S V' that is U diag(1, 1, d) V',
where d = det(U) det(V). The directions for the two largest singular values
are found here from the eigenvectors of k' k, and the third as their cross
product, on both sides. So the result is a proper rotation to within
rounding however small the third singular value is: reference directions
that all lie in one plane determine the rotation perfectly well.

Returns 0 on success, or -1 if `k` has (numerically) fewer than two
independent directions.
*/
int _trical_polar_rotation(const trical_wide_t k[9], trical_wide_t q[9]) {
    trical_wide_t ktk[9], w[3], v[9], u1[3], u2[3], u3[3], v1[3], v2[3];
    trical_wide_t v3[3], norm1, norm2, dot;
    unsigned int i, j, first = 0, second;

    assert(k && q);

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            ktk[i * 3 + j] = k[0 * 3 + i] * k[0 * 3 + j] +
                             k[1 * 3 + i] * k[1 * 3 + j] +
                             k[2 * 3 + i] * k[2 * 3 + j];
        }
    }

    trical_eigen_sym3(ktk, w, v);

    /* the two largest squared singular values */
    for (i = 1; i < 3; i++) {
        if (w[i] > w[first]) {
            first = i;
        }
    }
    second = (first + 1) % 3;
    if (w[(first + 2) % 3] > w[second]) {
        second = (first + 2) % 3;
    }

    for (i = 0; i < 3; i++) {
        v1[i] = v[i * 3 + first];
        v2[i] = v[i * 3 + second];
    }
    for (i = 0; i < 3; i++) {
        u1[i] = k[i * 3 + 0] * v1[0] + k[i * 3 + 1] * v1[1] +
                k[i * 3 + 2] * v1[2];
        u2[i] = k[i * 3 + 0] * v2[0] + k[i * 3 + 1] * v2[1] +
                k[i * 3 + 2] * v2[2];
    }

    norm1 = trical_sqrt(u1[0] * u1[0] + u1[1] * u1[1] + u1[2] * u1[2]);
    if (!(norm1 > 0) || !_wide_isfinite(norm1)) {
        return -1;
    }
    for (i = 0; i < 3; i++) {
        u1[i] /= norm1;
    }

    /* k v2 is orthogonal to k v1 already, but for rounding */
    dot = u1[0] * u2[0] + u1[1] * u2[1] + u1[2] * u2[2];
    for (i = 0; i < 3; i++) {
        u2[i] -= dot * u1[i];
    }
    norm2 = trical_sqrt(u2[0] * u2[0] + u2[1] * u2[1] + u2[2] * u2[2]);
    if (!(norm2 > norm1 * TRICAL_WIDE_EPSILON * 64) ||
            !_wide_isfinite(norm2)) {
        return -1;
    }
    for (i = 0; i < 3; i++) {
        u2[i] /= norm2;
    }

    u3[0] = u1[1] * u2[2] - u1[2] * u2[1];
    u3[1] = u1[2] * u2[0] - u1[0] * u2[2];
    u3[2] = u1[0] * u2[1] - u1[1] * u2[0];
    v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
    v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
    v3[2] = v1[0] * v2[1] - v1[1] * v2[0];

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            q[i * 3 + j] = u1[i] * v1[j] + u2[i] * v2[j] + u3[i] * v3[j];
        }
    }

    return 0;
}

/*
_trical_recentre
Moves the origin of the normalised coordinates by (approximately) `shift`,
re-expressing everything that is stored in those coordinates.
*/
void _trical_recentre(TRICAL_instance_t *instance,
const trical_wide_t shift[3]) {
    trical_wide_t delta[3];
    unsigned int i, j;

    for (i = 0; i < 3; i++) {
        /*
        The origin is stored in single precision, so apply exactly the shift
        that the rounded origin represents.
        */
        const float origin = (float)((trical_wide_t)instance->priv.origin[i] +
            (trical_wide_t)instance->priv.scale * shift[i]);

        delta[i] = ((trical_wide_t)origin -
                    (trical_wide_t)instance->priv.origin[i]) /
                   (trical_wide_t)instance->priv.scale;
        instance->priv.origin[i] = origin;
    }

    _trical_moments_shift(instance->priv.moments, delta);

    /* sum f (x - delta)' = sum f x' - (sum f) delta' */
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            instance->priv.cross[i * 3 + j] -=
                instance->priv.cross[9 + i] * delta[j];
        }
    }

    _trical_quadric_translate(instance->priv.quadric, delta);
    _trical_quadric_translate(instance->priv.quadric_prior, delta);

    for (i = 0; i < 3; i++) {
        instance->priv.data_mean[i] -= (float)delta[i];
    }
}

/*
Measurement noise standard deviation assumed by the fit, in normalised units
(i.e. as a fraction of the size of the ellipsoid): the estimate from the data
once there is one; until then the user's value, within limits.
*/
static trical_wide_t _noise_assumed(const TRICAL_instance_t *instance) {
    const trical_wide_t lower = (trical_wide_t)TRICAL_NOISE_INITIAL_MIN;
    const trical_wide_t upper = (trical_wide_t)TRICAL_NOISE_INITIAL_MAX;
    trical_wide_t noise;

    if (instance->priv.noise_estimate > 0) {
        return (trical_wide_t)instance->priv.noise_estimate;
    }

    noise = (trical_wide_t)instance->measurement_noise /
            (trical_wide_t)instance->priv.scale;

    /* written so that a NaN ends up at the lower limit */
    return noise > upper ? upper : (noise > lower ? noise : lower);
}

/*
Sets the size of the cells that measurements are grouped into (in the units
of the measurements): TRICAL_GATE times the size of the ellipsoid, or
TRICAL_GATE_NOISE times the noise if that is more, so that noise alone does
not carry the measurements of a sensor at rest from one cell to the next.
The number of measurements used per cell is in proportion to its size, so
that a moving sensor contributes the same number per distance moved.
*/
static void _gate_update(TRICAL_instance_t *instance) {
    const trical_wide_t noise = (trical_wide_t)TRICAL_GATE_NOISE *
                                _noise_assumed(instance);
    const trical_wide_t gate = noise > (trical_wide_t)TRICAL_GATE ? noise :
                               (trical_wide_t)TRICAL_GATE;
    const trical_wide_t quota = (trical_wide_t)TRICAL_CELL_QUOTA * gate /
                                (trical_wide_t)TRICAL_GATE;

    instance->priv.gate = (float)gate * instance->priv.scale;
    instance->priv.cell_quota = (unsigned short)(quota < 60000 ? quota :
                                                 60000);
}

/*
_trical_start
Starts accumulating, taking the current `state` of `instance` as the initial
estimate -- except that the initial bias is `bias`, if that is not NULL.
`state` itself is not changed. Returns 0 on success, or -1 (leaving
`instance` untouched) if `state` is not a usable calibration.
*/
int _trical_start(TRICAL_instance_t *instance, const float bias[3]) {
    const trical_wide_t ratio2 = (trical_wide_t)(TRICAL_AXIS_RATIO_MAX *
                                                 TRICAL_AXIS_RATIO_MAX);
    trical_wide_t a[9], ata[9], w[3], v[9], d[3], s[9], rotation[9];
    trical_wide_t centre[3] = { 0, 0, 0 }, wmin, wmax, mean_scale = 0;
    float origin[3];
    unsigned int i, j;

    assert(instance);

    /* (I + D) from the state; reject anything that is not finite */
    for (i = 0; i < 3; i++) {
        origin[i] = bias ? bias[i] : instance->state[i];
        if (!_wide_isfinite((trical_wide_t)origin[i])) {
            return -1;
        }
        for (j = 0; j < 3; j++) {
            a[i * 3 + j] = (trical_wide_t)instance->state[3 + i * 3 + j] +
                           (i == j ? 1 : 0);
            if (!_wide_isfinite(a[i * 3 + j])) {
                return -1;
            }
        }
    }
    if (!_float_isfinite(instance->field_norm) ||
            !(instance->field_norm >= TRICAL_NORM_MIN) ||
            !(instance->field_norm <= TRICAL_NORM_MAX)) {
        return -1;
    }

    /* (I + D)' (I + D) defines the ellipsoid; it must be a sane one */
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            ata[i * 3 + j] = a[0 * 3 + i] * a[0 * 3 + j] +
                             a[1 * 3 + i] * a[1 * 3 + j] +
                             a[2 * 3 + i] * a[2 * 3 + j];
        }
    }
    trical_eigen_sym3(ata, w, v);

    wmin = w[0] < w[1] ? (w[0] < w[2] ? w[0] : w[2]) :
                         (w[1] < w[2] ? w[1] : w[2]);
    wmax = w[0] > w[1] ? (w[0] > w[2] ? w[0] : w[2]) :
                         (w[1] > w[2] ? w[1] : w[2]);
    if (!(wmin > 0) || !(wmax <= ratio2 * wmin) || !_wide_isfinite(wmax)) {
        return -1;
    }

    /* rotation part: (I + D) [(I + D)' (I + D)]^(-1/2) */
    for (i = 0; i < 3; i++) {
        d[i] = (trical_wide_t)1 / trical_sqrt(w[i]);
        mean_scale += trical_sqrt(w[i]) / 3;
    }
    _similarity3(v, d, s);
    _mul33(a, s, rotation);

    /*
    If that is the identity to within rounding (as it is for any symmetric
    initial estimate), make it exactly so: the scale estimate is then
    exactly symmetric until reference directions say otherwise.
    */
    wmax = 0;
    for (i = 0; i < 9; i++) {
        const trical_wide_t difference =
            trical_abs(rotation[i] - (i % 4 == 0 ? 1 : 0));

        if (difference > wmax) {
            wmax = difference;
        }
    }
    if (wmax < (trical_wide_t)1e-5) {
        for (i = 0; i < 9; i++) {
            rotation[i] = (i % 4 == 0) ? 1 : 0;
        }
    }

    /*
    The size of the ellipsoid in the units of the measurements, which is
    the unit of the normalised coordinates, has to be such that squares of
    distances between measurements are representable.
    */
    wmin = (trical_wide_t)instance->field_norm / mean_scale;
    if (!(wmin >= (trical_wide_t)TRICAL_SCALE_MIN) ||
            !(wmin <= (trical_wide_t)TRICAL_SCALE_MAX)) {
        return -1;
    }

    /* all checks passed: commit */
    memset(&instance->priv, 0, sizeof(instance->priv));

    for (i = 0; i < 9; i++) {
        instance->priv.rotation_prior[i] = (float)rotation[i];
    }
    memcpy(instance->priv.origin, origin, sizeof(instance->priv.origin));
    instance->priv.data_range2 = TRICAL_RANGE * TRICAL_RANGE;

    /*
    Normalised coordinates: x = (measurement - origin) / scale, with the
    origin at the initial bias and the scale such that the initial ellipsoid
    |(I + D) scale x| = field norm has a mean radius of about 1.
    */
    instance->priv.scale = (float)wmin;
    for (i = 0; i < 9; i++) {
        ata[i] /= mean_scale * mean_scale;
    }
    _trical_quadric_from_calibration(ata, centre,
                                     instance->priv.quadric_prior);
    memcpy(instance->priv.quadric, instance->priv.quadric_prior,
           sizeof(instance->priv.quadric));

    memcpy(instance->priv.published, instance->state,
           sizeof(instance->priv.published));
    instance->priv.started = 1;
    _gate_update(instance);

    return 0;
}

/*
Centre c = -G^-1 q, squared size rho = q' G^-1 q - q0 and trace of G of the
ellipsoid `theta`, with the sign of theta taken such that the trace is
positive. Returns 0 on success, or -1 if `theta` is not an ellipsoid.
*/
static int _ellipsoid(const trical_wide_t theta[10], trical_wide_t centre[3],
trical_wide_t *rho, trical_wide_t *trace) {
    trical_wide_t g[9], q[3], sign = 1;
    unsigned int i;

    *trace = theta[0] + theta[1] + theta[2];
    if (*trace < 0) {
        sign = -1;
        *trace = -*trace;
    }

    /* lower triangle of G */
    g[0] = sign * theta[0];
    g[3] = sign * theta[5];
    g[4] = sign * theta[1];
    g[6] = sign * theta[4];
    g[7] = sign * theta[3];
    g[8] = sign * theta[2];
    if (trical_cholesky(g, 3) != 0) {
        return -1;
    }

    for (i = 0; i < 3; i++) {
        q[i] = sign * theta[6 + i];
        centre[i] = -q[i];
    }
    trical_cholesky_solve(g, 3, centre);

    *rho = -(q[0] * centre[0] + q[1] * centre[1] + q[2] * centre[2]) -
           sign * theta[9];
    if (!(*rho > 0) || !_wide_isfinite(*rho)) {
        return -1;
    }

    return 0;
}

/*
The prior, which holds the fit to the initial estimate where the
measurements do not determine it. It is a quadratic penalty theta' Q theta,
with Q the sum of three positive semi-definite terms:

- shape (rank 5): the part of G that is not proportional to the G of the
  prior quadric, i.e. differences in the ratios and directions of the
  principal scale factors;
- size (rank 1): (tr G - k rho)^2 / 3, where k is tr G / rho of the prior
  quadric, i.e. a difference in the overall scale factor. rho is not linear
  in theta, so -F(c) is used instead, with c the centre of the current
  solution, which agrees with it to second order in the distance of c from
  the centre of the solution;
- bias (rank 3): |G c0 + q|^2 = |G (c0 - c)|^2, where c0 is the centre of
  the prior quadric and c that of the solution.

With G = t (I + 2 dA) for small differences dA in the scale factors, these
are 4 t^2 |dA|^2 (Frobenius) and t^2 |c - c0|^2; and the squared distance of
a measurement from the ellipsoid, d^2, contributes 4 t^2 d^2 to
theta' M theta. So with a weight of w the prior is worth w (|dA|^2 +
|c - c0|^2 / 4) of mean squared distance, whatever the coordinates.

The first two terms are what makes sense of measurements that cover only
part of the ellipsoid: the scale factors are seldom far from their initial
values. The bias can be, by several times the size of the ellipsoid, and a
quadratic penalty on that would then swamp everything. So the bias term is
weighted down by 1 / (1 + d^2 / TRICAL_PRIOR_BIAS_RANGE^2), where d is the
RMS distance of the measurements from the prior ellipsoid: the evidence
that the initial estimate is wrong. (That depends only on the measurements,
not on the current solution, so noise cannot wear the hold on the bias down
by pushing the solution about.)

The shape and size terms can be strengthened by a factor, which in the limit
gives the fit of an ellipsoid of the initial shape and size wherever the
measurements put it -- which is always an acceptable ellipsoid.
*/
typedef struct {
    /* shape: the unit vector, in the Frobenius metric, of the prior G */
    trical_wide_t shape[6];
    /* size and (weighted) bias terms, as row vectors */
    trical_wide_t size[10];
    trical_wide_t bias[3][10];
    /* factors by which the shape and size terms are strengthened */
    trical_wide_t shape_boost;
    trical_wide_t size_boost;
    /* tr G / rho of the prior quadric */
    trical_wide_t kappa;
    /*
    how well the prior quadric itself fits the measurements: 1 if as well
    as the current solution does, falling to 0 as it is found wanting
    */
    trical_wide_t hold;
} _prior_t;

/* weights of the elements of theta in the Frobenius norm of G */
static const unsigned char TRICAL_G_METRIC[6] = { 1, 1, 1, 2, 2, 2 };

/*
Sets the size term of `prior` for a solution whose centre is (near) `c`:
(tr G + kappa F(c)) / sqrt(3), with F(c) = xi(c) . theta.
*/
static void _prior_size(_prior_t *prior, const trical_wide_t c[3]) {
    const trical_wide_t root = trical_sqrt((trical_wide_t)3);
    const trical_wide_t kappa = prior->kappa;

    prior->size[0] = (kappa * c[0] * c[0] + 1) / root;
    prior->size[1] = (kappa * c[1] * c[1] + 1) / root;
    prior->size[2] = (kappa * c[2] * c[2] + 1) / root;
    prior->size[3] = kappa * 2 * c[1] * c[2] / root;
    prior->size[4] = kappa * 2 * c[0] * c[2] / root;
    prior->size[5] = kappa * 2 * c[0] * c[1] / root;
    prior->size[6] = kappa * 2 * c[0] / root;
    prior->size[7] = kappa * 2 * c[1] / root;
    prior->size[8] = kappa * 2 * c[2] / root;
    prior->size[9] = kappa / root;
}

/*
Sets up the prior of `instance` about its current solution, given the mean
quadric basis `xi` of the measurements. Returns 0 on success, or -1 if the
prior or current quadric is not an ellipsoid (which cannot happen unless
something is corrupt).
*/
static int _prior_init(const TRICAL_instance_t *instance,
const trical_wide_t xi[10], _prior_t *prior) {
    const trical_wide_t *theta0 = instance->priv.quadric_prior;
    const trical_wide_t range2 = (trical_wide_t)(TRICAL_PRIOR_BIAS_RANGE *
                                                 TRICAL_PRIOR_BIAS_RANGE);
    const trical_wide_t count = instance->priv.moments[0];
    trical_wide_t c0[3], c[3], rho0, rho, trace0, trace, norm = 0;
    trical_wide_t mt[10], nt[10], fit, gradient2;
    trical_wide_t kappa, beta;
    unsigned int i;

    if (_ellipsoid(theta0, c0, &rho0, &trace0) != 0 ||
            _ellipsoid(instance->priv.quadric, c, &rho, &trace) != 0) {
        return -1;
    }
    prior->shape_boost = 1;
    prior->size_boost = 1;

    for (i = 0; i < 6; i++) {
        norm += (trical_wide_t)TRICAL_G_METRIC[i] * theta0[i] * theta0[i];
    }
    norm = trical_sqrt(norm);
    for (i = 0; i < 6; i++) {
        prior->shape[i] = theta0[i] / norm;
    }

    kappa = trace0 / rho0;
    prior->kappa = kappa;
    _prior_size(prior, c);

    /*
    The evidence that the initial bias is wrong sets the weight of the bias
    term: the mean squared distance of the measurements from the prior
    ellipsoid (Taubin's approximation, which reads low far from the
    ellipsoid); or the square of the distance by which the mean of the
    measurements lies outside it, if that is more (no measurement can be
    further from the bias than the size of the ellipsoid, 3 / kappa being
    its mean squared semi-axis). Written so that a NaN, or measurements
    where Taubin's approximation breaks down (at the centre of the
    ellipsoid), give no weight.
    */
    _trical_moment_product(instance->priv.moments, theta0, mt);
    _trical_normalisation_product(xi, theta0, 0, nt);
    fit = _dot10(theta0, mt);
    gradient2 = _dot10(theta0, nt);

    /* rounding can leave a perfect fit slightly negative */
    if (fit < 0) {
        fit = 0;
    }

    beta = 0;
    prior->hold = 0;
    if (gradient2 > 0 && fit >= 0) {
        trical_wide_t evidence2 = fit / gradient2, excess = 0;

        for (i = 0; i < 3; i++) {
            const trical_wide_t mean = xi[6 + i] / 2 - c0[i];

            excess += mean * mean;
        }
        excess = trical_sqrt(excess) - trical_sqrt(3 / kappa);
        if (excess > 0 && excess * excess > evidence2) {
            evidence2 = excess * excess;
        }

        beta = range2 / (range2 + evidence2);
        beta = (beta >= 0) ? trical_sqrt(beta) : 0;

        /*
        Is there, on the other hand, any evidence against the initial
        estimate at all? Not if the prior ellipsoid fits the measurements
        as well as the current solution does, which has had nine parameters
        adjusted to them: then whatever separates the two is something
        that the measurements do not determine, and the prior should hold
        it firmly (see _regularisation).
        */
        if (count > 9) {
            const trical_wide_t limit =
                (trical_wide_t)TRICAL_PRIOR_HOLD_EXCESS;
            trical_wide_t best, floor2;

            /*
            The mean squared distance of the measurements from an ellipsoid
            cannot be told from zero below the noise floor, or below the
            rounding errors of the sums.
            */
            floor2 = (trical_wide_t)TRICAL_PRIOR_HOLD_ROUNDING *
                     trical_sqrt(count) * TRICAL_WIDE_EPSILON;
            if (floor2 < (trical_wide_t)(TRICAL_NOISE_FLOOR *
                                         TRICAL_NOISE_FLOOR)) {
                floor2 = (trical_wide_t)(TRICAL_NOISE_FLOOR *
                                         TRICAL_NOISE_FLOOR);
            }

            _trical_moment_product(instance->priv.moments,
                                   instance->priv.quadric, mt);
            _trical_normalisation_product(xi, instance->priv.quadric, 0, nt);
            fit = _dot10(instance->priv.quadric, mt);
            gradient2 = _dot10(instance->priv.quadric, nt);

            if (gradient2 > 0) {
                best = (fit > 0 ? fit / gradient2 : 0) * count / (count - 9);
                if (!(best > floor2)) {
                    best = floor2;
                }

                excess = evidence2 / best - 1;
                if (excess > 0) {
                    excess /= limit;
                    prior->hold = 1 / (1 + excess * excess);
                } else if (excess <= 0) {
                    prior->hold = 1;
                }
            }
        }
    }

    /* sqrt(beta) (G c0 + q) */
    memset(prior->bias, 0, sizeof(prior->bias));
    prior->bias[0][0] = beta * c0[0];
    prior->bias[0][5] = beta * c0[1];
    prior->bias[0][4] = beta * c0[2];
    prior->bias[0][6] = beta;
    prior->bias[1][5] = beta * c0[0];
    prior->bias[1][1] = beta * c0[1];
    prior->bias[1][3] = beta * c0[2];
    prior->bias[1][7] = beta;
    prior->bias[2][4] = beta * c0[0];
    prior->bias[2][3] = beta * c0[1];
    prior->bias[2][2] = beta * c0[2];
    prior->bias[2][8] = beta;

    return 0;
}

/* y = Q theta for the prior matrix Q */
static void _prior_product(const _prior_t *prior,
const trical_wide_t theta[10], trical_wide_t y[10]) {
    trical_wide_t projection = 0, size = _dot10(prior->size, theta);
    unsigned int i, j;

    for (i = 0; i < 6; i++) {
        projection += (trical_wide_t)TRICAL_G_METRIC[i] * prior->shape[i] *
                      theta[i];
    }
    for (i = 0; i < 10; i++) {
        y[i] = prior->size_boost * size * prior->size[i];
        if (i < 6) {
            y[i] += prior->shape_boost * (trical_wide_t)TRICAL_G_METRIC[i] *
                    (theta[i] - projection * prior->shape[i]);
        }
    }
    for (j = 0; j < 3; j++) {
        const trical_wide_t bias = _dot10(prior->bias[j], theta);

        for (i = 0; i < 10; i++) {
            y[i] += bias * prior->bias[j][i];
        }
    }
}

/* Element (i, j) of the prior matrix Q */
static trical_wide_t _prior_element(const _prior_t *prior, unsigned int i,
unsigned int j) {
    trical_wide_t value = prior->size_boost * prior->size[i] * prior->size[j] +
                          prior->bias[0][i] * prior->bias[0][j] +
                          prior->bias[1][i] * prior->bias[1][j] +
                          prior->bias[2][i] * prior->bias[2][j];

    if (i < 6 && j < 6) {
        trical_wide_t shape = -(trical_wide_t)(TRICAL_G_METRIC[i] *
            TRICAL_G_METRIC[j]) * prior->shape[i] * prior->shape[j];

        if (i == j) {
            shape += (trical_wide_t)TRICAL_G_METRIC[i];
        }
        value += prior->shape_boost * shape;
    }

    return value;
}

/*
The Rayleigh quotient theta' M theta / theta' N theta: an estimate of the
variance of the measurement noise, if theta fits the measurements as well as
an ellipsoid can. Zero if it is not positive.
*/
static trical_wide_t _quotient(const trical_wide_t *sums,
const trical_wide_t xi[10], const trical_wide_t theta[10]) {
    trical_wide_t mt[10], nt[10], fit, norm;

    _trical_moment_product(sums, theta, mt);
    _trical_normalisation_product(xi, theta, 1, nt);
    fit = _dot10(theta, mt);
    norm = _dot10(theta, nt);

    return (fit > 0 && norm > 0) ? fit / norm : 0;
}

/*
Works out, from the current solution, the two numbers that define the next
eigenproblem (see _trical_refit): the noise correction `lambda` and the
weight of the prior.
*/
static void _regularisation(const TRICAL_instance_t *instance,
const trical_wide_t xi[10], const _prior_t *prior, trical_wide_t *lambda,
trical_wide_t *weight) {
    const trical_wide_t *theta = instance->priv.quadric;
    const trical_wide_t count = instance->priv.moments[0];
    const trical_wide_t noise = _noise_assumed(instance);
    trical_wide_t limit;

    /*
    lambda: the Rayleigh quotient of the current solution, which converges
    on the smallest eigenvalue of M theta = lambda N theta. It is an estimate
    of the noise variance, so it is limited to what that could credibly be:
    until there are enough measurements to tell, a few times what the noise
    is assumed to be (the first few measurements may be nowhere near the
    initial estimate, and that is not noise); after that, a tenth of the
    size of the ellipsoid -- correcting for more would not rescue a fit that
    poor. Written so that a NaN ends up at the limit.
    */
    limit = (trical_wide_t)TRICAL_CORRECTION_MAX;
    if (count < (trical_wide_t)TRICAL_NOISE_MIN_COUNT &&
            4 * noise * noise < limit) {
        limit = 4 * noise * noise;
    }

    *lambda = _quotient(instance->priv.moments, xi, theta);
    if (!(*lambda < limit)) {
        *lambda = limit;
    }

    /*
    The weight of the prior, in units of squared distance:
        noise^2 (TRICAL_PRIOR_BASE + g sqrt(count)) / count
    The first term is a prior standard deviation of about a tenth in the
    scale factors, and a tenth of the field norm in the bias, for
    TRICAL_PRIOR_BASE = 300.

    The second is there because noise does not leave alone what the
    measurements do not determine: the fit of an ellipsoid to noisy
    measurements can always be improved a little, by chance, by moving it
    in such a direction, and the improvement on offer grows as sqrt(count).
    The prior has to outweigh that, or the estimate wanders -- by about
    4 / g, as a fraction of the field norm. So g is TRICAL_PRIOR_GROWTH,
    which is as much as an initial estimate that the measurements contradict
    can be allowed without getting in their way, plus TRICAL_PRIOR_HOLD to
    the extent that they do not contradict it.
    */
    *weight = noise * noise *
        ((trical_wide_t)TRICAL_PRIOR_BASE +
         ((trical_wide_t)TRICAL_PRIOR_GROWTH +
          (trical_wide_t)TRICAL_PRIOR_HOLD * prior->hold) *
         trical_sqrt(count)) / count;

    /*
    ...but never less than the rounding errors of the sums amount to: with
    a margin, to the extent that the initial estimate is worth holding on
    to.
    */
    limit = (1 + ((trical_wide_t)TRICAL_PRIOR_ROUNDING - 1) * prior->hold) *
            trical_sqrt(count) * TRICAL_WIDE_EPSILON;
    if (*weight < limit) {
        *weight = limit;
    }
}

/*
Writes to `c` the Cholesky factor of
    M - lambda N + weight Q + shift I
where M and N are the moment and (hyper) normalisation matrices of the
measurements and Q is the prior matrix, where the shift makes the matrix
safely positive definite: `shift` to begin with, and ten times more at each
further attempt, up to TRICAL_SHIFT_MAX; the shift that worked is left in
`shift`. Each factorisation attempted takes one from `budget`. Returns 0 on
success, or -1 if no reasonable shift does or the budget runs out.
*/
static int _factorise(const trical_wide_t *sums, const trical_wide_t xi[10],
const _prior_t *prior, trical_wide_t lambda, trical_wide_t weight,
trical_wide_t c[100], trical_wide_t *shift, unsigned int *budget) {
    const trical_wide_t scale = (trical_wide_t)1 / sums[0];
    unsigned int i, j, k;

    if (!(*shift > TRICAL_SHIFT)) {
        *shift = TRICAL_SHIFT;
    }

    for (; *shift < (trical_wide_t)TRICAL_SHIFT_MAX && *budget > 0;
            *shift *= 10) {
        *budget -= 1;

        /* only the lower triangle is needed */
        for (i = 0; i < 10; i++) {
            for (j = 0; j <= i; j++) {
                c[i * 10 + j] =
                    (trical_wide_t)(TRICAL_XI_COEF[i] * TRICAL_XI_COEF[j]) *
                    sums[TRICAL_M_INDEX[i][j]] * scale +
                    weight * _prior_element(prior, i, j);
            }
            c[i * 10 + i] += *shift;
        }

        /* - lambda N_T (its terms are of the upper triangle) */
        for (k = 0; k < 30; k++) {
            c[TRICAL_N_TERMS[k][1] * 10 + TRICAL_N_TERMS[k][0]] -= lambda *
                (trical_wide_t)TRICAL_N_TERMS[k][3] * xi[TRICAL_N_TERMS[k][2]];
        }
        /* - lambda (xi e' + e xi') */
        for (j = 0; j < 3; j++) {
            for (i = j; i < 10; i++) {
                c[i * 10 + j] -= lambda * (xi[i] + (i < 3 ? xi[j] : 0));
            }
        }

        if (trical_cholesky(c, 10) == 0) {
            return 0;
        }
    }

    return -1;
}

/*
Finds the eigenvector for the smallest eigenvalue of the matrix whose
(shifted) Cholesky factor is `c`, by inverse iteration from `theta`. Returns
the squared change in theta over the last iteration, or -1 if the iteration
broke down.
*/
static trical_wide_t _smallest_eigenvector(const trical_wide_t c[100],
trical_wide_t theta[10]) {
    /* squared change in theta at which the iteration has converged */
    const trical_wide_t tolerance =
        (sizeof(trical_wide_t) == sizeof(float)) ? (trical_wide_t)1e-8 :
                                                   (trical_wide_t)1e-14;
    trical_wide_t y[10], change = 0;
    unsigned int i, iteration;

    for (iteration = 0; iteration < TRICAL_ITERATIONS_MAX; iteration++) {
        trical_wide_t norm;

        memcpy(y, theta, sizeof(y));
        trical_cholesky_solve(c, 10, y);

        /* theta' C^-1 theta > 0, so the sign of theta is preserved */
        norm = trical_sqrt(_dot10(y, y));
        if (!(norm > 0) || !_wide_isfinite(norm)) {
            return -1;
        }

        change = 0;
        for (i = 0; i < 10; i++) {
            const trical_wide_t next = y[i] / norm;

            change += (next - theta[i]) * (next - theta[i]);
            theta[i] = next;
        }

        if (iteration + 1 >= TRICAL_ITERATIONS_MIN && change < tolerance) {
            break;
        }
    }

    return change;
}

/*
Writes the calibration (centre, a) of the normalised ellipsoid to the state
of `instance`, including the rotation derived from the reference directions
-- but leaving the bias as it is if `keep_bias` is non-zero. Returns 0 if it
did, or -1 if the result cannot be represented (which takes an absurd field
norm), in which case the state is left as it was.
*/
static int _publish(TRICAL_instance_t *instance,
const trical_wide_t centre[3], const trical_wide_t a[9], int keep_bias) {
    const trical_wide_t gain = (trical_wide_t)instance->field_norm /
                               (trical_wide_t)instance->priv.scale;
    trical_wide_t rotation[9], k[9], full[9];
    float next[TRICAL_STATE_DIM];
    unsigned int i, j;
    int aligned = -1;

    if (instance->priv.reference_count > 0) {
        /*
        Attitude profile matrix K = E[f c'], where c = a (x - centre) is the
        calibrated reading scaled to unit length:
            sum f c' = (sum f x' - (sum f) centre') a
        A small multiple of the prior rotation settles whatever the data
        leave undetermined.
        */
        const trical_wide_t *cross = instance->priv.cross;
        const trical_wide_t scale = (trical_wide_t)1 /
            (trical_wide_t)instance->priv.reference_count;
        trical_wide_t fx[9];

        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                fx[i * 3 + j] = (cross[i * 3 + j] -
                                 cross[9 + i] * centre[j]) * scale;
            }
        }
        _mul33(fx, a, k);
        for (i = 0; i < 9; i++) {
            k[i] += (trical_wide_t)TRICAL_ROTATION_PRIOR *
                    (trical_wide_t)instance->priv.rotation_prior[i];
        }

        aligned = _trical_polar_rotation(k, rotation);
    }
    if (aligned != 0) {
        for (i = 0; i < 9; i++) {
            rotation[i] = (trical_wide_t)instance->priv.rotation_prior[i];
        }
    }

    _mul33(rotation, a, full);

    /* whatever else happens, nothing that is not a number gets out */
    for (i = 0; i < 3; i++) {
        trical_wide_t value = (trical_wide_t)instance->priv.origin[i] +
                              (trical_wide_t)instance->priv.scale * centre[i];

        if (!_wide_isfloat(value)) {
            return -1;
        }
        next[i] = keep_bias ? instance->state[i] : (float)value;

        for (j = 0; j < 3; j++) {
            value = gain * full[i * 3 + j] - (i == j ? 1 : 0);
            if (!_wide_isfloat(value)) {
                return -1;
            }
            next[3 + i * 3 + j] = (float)value;
        }
    }

    memcpy(instance->state, next, sizeof(instance->state));
    memcpy(instance->priv.published, next, sizeof(instance->priv.published));

    return 0;
}

/*
One step of the fit: finds the unit vector theta that makes

    theta' (M - lambda N + w Q) theta

smallest, where w Q is the prior about the current solution of `instance`
(see _prior_t), and lambda is `correction`, or the Rayleigh quotient
theta' M theta / theta' N theta of the current solution if `correction` is
negative (in either case subject to the limits in _regularisation). That is
an ordinary symmetric eigenproblem, solved by inverse iteration from the
current solution.

`budget` is the number of matrix factorisations (the bulk of the work) that
may still be made for this refit; the step fails rather than exceed it. That
bounds the time a refit can take, whatever the data.

The size term of the prior presumes that the centre of the result is near
that of the current solution. If it turns out not to be, the term is set up
again about the new centre and the problem solved again, a few times if
need be. If the result is not an ellipsoid at all, the centre of the best
fit of an ellipsoid of the initial shape (which is a much better-posed
problem) is tried instead.

If the result is still not an acceptable ellipsoid, the step falls back on
the one without the noise correction, which is biased but as stable as can
be; and if that will not do either, the hold on the scale factors is
strengthened step by step, which leads to an ellipsoid of the initial shape
and size fitted to the measurements by its position alone.

Returns 0 on success, with the solution in `theta` and as (`centre`, `a`),
the Cholesky factor of the matrix in `c`, the correction that was used in
`lambda`, and the squared change in theta over the last iteration in
`change`; or -1 on failure.
*/
static int _fit_step(const TRICAL_instance_t *instance,
const trical_wide_t xi[10], trical_wide_t correction, _prior_t *prior,
trical_wide_t c[100],
trical_wide_t theta[10], trical_wide_t centre[3], trical_wide_t a[9],
trical_wide_t *lambda, trical_wide_t *weight, trical_wide_t *change,
trical_wide_t *shift, unsigned int *budget) {
    const trical_wide_t limit2 = (trical_wide_t)(TRICAL_PRIOR_RECENTRE *
                                                 TRICAL_PRIOR_RECENTRE);
    const trical_wide_t *sums = instance->priv.moments;
    trical_wide_t base[3], about[3], found[3], rho, trace, probe_shift = 0;
    unsigned int i, attempt = 0, recentred = 0;
    int probed = 0;

    if (_prior_init(instance, xi, prior) != 0 ||
            _ellipsoid(instance->priv.quadric, base, &rho, &trace) != 0) {
        return -1;
    }
    memcpy(about, base, sizeof(about));
    _regularisation(instance, xi, prior, lambda, weight);
    if (correction >= 0 && correction < *lambda) {
        *lambda = correction;
    }

    for (;;) {
        int ellipsoid = 0;

        if (*budget == 0) {
            return -1;
        }

        *change = -1;
        if (_factorise(sums, xi, prior, *lambda, *weight, c, shift,
                       budget) == 0) {
            memcpy(theta, instance->priv.quadric, 10 * sizeof(theta[0]));
            *change = _smallest_eigenvector(c, theta);
        }
        if (*change >= 0) {
            ellipsoid = _ellipsoid(theta, found, &rho, &trace) == 0;
        }

        if (!ellipsoid && !probed) {
            /*
            Where would an ellipsoid of the initial shape, of any size, be?
            (Strengthening the hold on the shape to the exclusion of
            everything else gives that.)
            */
            const trical_wide_t shape = prior->shape_boost,
                                size = prior->size_boost;

            probed = 1;
            prior->shape_boost = (trical_wide_t)TRICAL_FIT_PROBE / *weight;
            prior->size_boost = 0;
            if (_factorise(sums, xi, prior, 0, *weight, c, &probe_shift,
                           budget) == 0) {
                memcpy(theta, instance->priv.quadric,
                       10 * sizeof(theta[0]));

                /*
                Only an ellipsoid of a plausible size says anything about
                where the centre is: measurements that are all in one
                place, for instance, are best fitted by one that has shrunk
                to nothing around them.
                */
                if (_smallest_eigenvector(c, theta) >= 0 &&
                        _trical_quadric_to_calibration(theta, found, a,
                                                       NULL) == 0) {
                    ellipsoid = 1;
                }
            }
            prior->shape_boost = shape;
            prior->size_boost = size;

            if (ellipsoid) {
                memcpy(base, found, sizeof(base));
                memcpy(about, base, sizeof(about));
                _prior_size(prior, about);
                continue;
            }
        }

        /* is the size term set up about the right place? */
        if (ellipsoid && recentred < TRICAL_PRIOR_RECENTRE_LIMIT) {
            trical_wide_t distance2 = 0;

            for (i = 0; i < 3; i++) {
                distance2 += (found[i] - about[i]) * (found[i] - about[i]);
            }
            if (!(distance2 < limit2)) {
                memcpy(about, found, sizeof(about));
                _prior_size(prior, about);
                recentred++;
                continue;
            }
        }

        if (ellipsoid &&
                _trical_quadric_to_calibration(theta, centre, a, NULL) == 0) {
            return 0;
        }
        if (attempt++ >= TRICAL_FIT_ATTEMPTS) {
            return -1;
        }
        if (*lambda > 0) {
            /* without the correction the matrix is positive semi-definite */
            *lambda = 0;
            *shift = 0;
        } else {
            prior->shape_boost *= (trical_wide_t)TRICAL_FIT_PRIOR_STEP;
            prior->size_boost *= (trical_wide_t)TRICAL_FIT_PRIOR_STEP;
        }

        /* start again from where the size term was first set up */
        memcpy(about, base, sizeof(about));
        _prior_size(prior, about);
    }
}

/*
_trical_republish
Writes to `state` the calibration that the current solution of `instance`
amounts to with its current field norm. Returns 0 on success, or -1 (and
`state` is untouched) if that is not a set of finite numbers.
*/
int _trical_republish(TRICAL_instance_t *instance) {
    trical_wide_t centre[3], a[9];

    assert(instance);

    if (!instance->priv.started ||
            _trical_quadric_to_calibration(instance->priv.quadric, centre, a,
                                           NULL) != 0) {
        return -1;
    }

    /*
    Until there has been a fit, the solution is the initial estimate, whose
    bias is the one in `state` -- or, after a restart, a first guess that
    is not worth publishing: `state` keeps the bias it has.
    */
    return _publish(instance, centre, a,
                    memcmp(instance->priv.quadric,
                           instance->priv.quadric_prior,
                           sizeof(instance->priv.quadric)) == 0);
}

/*
_trical_refit
Refits the calibration to the measurements accumulated so far, and updates
`state` if a valid fit is obtained.

Each refit is one step of the "renormalisation" iteration of [1] (see
_fit_step), which converges on the solution of M theta = lambda N theta as
fits succeed one another. That relies on the previous solution being a fair
approximation to the next one. Where it evidently is not -- the fit comes
out much closer to the measurements than the noise it was corrected for --
the step is repeated with the correction that the fit itself implies.
*/
void _trical_refit(TRICAL_instance_t *instance) {
    trical_wide_t c[100];
    trical_wide_t xi[10], theta[10], centre[3], a[9];
    trical_wide_t count, overlap, lambda, weight, change, shift;
    _prior_t prior;
    unsigned int i, budget = TRICAL_FIT_BUDGET;

    assert(instance);

    count = instance->priv.moments[0];
    if (!instance->priv.started || !(count > 0)) {
        return;
    }

    /*
    Sums in single precision cannot grow for ever: beyond a certain number
    of measurements, adding one more is lost in rounding. Halving the sums
    (and with them the weight of everything so far) before that happens
    costs nothing but a very long memory.
    */
    if ((trical_wide_t)TRICAL_COUNT_MAX > 0 &&
            count >= (trical_wide_t)TRICAL_COUNT_MAX) {
        for (i = 0; i < TRICAL_NUM_MOMENTS; i++) {
            instance->priv.moments[i] /= 2;
        }
        for (i = 0; i < 12; i++) {
            instance->priv.cross[i] /= 2;
        }
        instance->priv.reference_count -= instance->priv.reference_count / 2;
        count = instance->priv.moments[0];
    }

    /*
    Keep the origin of the normalised coordinates near the data, which keeps
    the problem well conditioned.
    */
    {
        trical_wide_t mean[3], distance2 = 0;

        for (i = 0; i < 3; i++) {
            mean[i] = instance->priv.moments[1 + i] / count;
            distance2 += mean[i] * mean[i];
        }
        if (distance2 > (trical_wide_t)(TRICAL_RECENTRE * TRICAL_RECENTRE)) {
            _trical_recentre(instance, mean);
        }
    }

    _trical_basis_mean(instance->priv.moments, xi);

    /*
    Refresh the range of measurements accepted as plausible: within a few
    times the RMS radius of the data from their mean (xi holds the mean
    squares of the coordinates, and twice their means).
    */
    if (count >= (trical_wide_t)TRICAL_OUTLIER_MIN_COUNT) {
        trical_wide_t radius2 = xi[0] + xi[1] + xi[2];
        float range;

        for (i = 0; i < 3; i++) {
            const trical_wide_t mean = xi[6 + i] / 2;

            radius2 -= mean * mean;
            instance->priv.data_mean[i] = (float)mean;
        }
        range = TRICAL_OUTLIER_RANGE *
                (1.0f + (float)trical_sqrt(radius2 > 0 ? radius2 : 0));
        instance->priv.data_range2 = range * range;
    }

    /* begin a little below the shift that the last fit needed */
    shift = (trical_wide_t)instance->priv.shift / 10;

    if (_fit_step(instance, xi, -1, &prior, c, theta, centre, a, &lambda,
                  &weight, &change, &shift, &budget) != 0) {
        return;
    }
    if (lambda > 0) {
        const trical_wide_t quotient = _quotient(instance->priv.moments, xi,
                                                 theta);

        if (quotient < (trical_wide_t)TRICAL_CORRECTION_REPEAT * lambda) {
            if (_fit_step(instance, xi, quotient, &prior, c, theta, centre,
                          a, &lambda, &weight, &change, &shift,
                          &budget) != 0) {
                return;
            }
        }
    }
    instance->priv.shift = (float)shift;

    overlap = trical_abs(_dot10(theta, instance->priv.quadric));
    memcpy(instance->priv.quadric, theta, sizeof(theta));

    /*
    Estimate the measurement noise from the fit: the Taubin quotient
    theta' M theta / theta' N_T theta is the mean squared distance of the
    measurements from the ellipsoid.

    If the fit is being held away from the data by the prior, or has not
    settled yet, that distance includes misfit as well as noise, so it is an
    upper bound. Lowering the assumed noise to it is therefore always safe.
    Raising the assumed noise is only done when the fit has settled and
    would not be much better with a much weaker prior -- otherwise a poor
    fit would be mistaken for noise, which would strengthen the prior, which
    would keep the fit poor.
    */
    if (count >= (trical_wide_t)TRICAL_NOISE_MIN_COUNT) {
        trical_wide_t mt[10], nt[10], fit, norm_t, estimate = -1;

        _trical_moment_product(instance->priv.moments, theta, mt);
        _trical_normalisation_product(xi, theta, 0, nt);
        fit = _dot10(theta, mt);
        norm_t = _dot10(theta, nt);

        /* rounding can leave a perfect fit slightly negative */
        if (!(fit > 0)) {
            fit = 0;
        }

        if (norm_t > 0) {
            /* nine parameters were fitted */
            estimate = trical_sqrt(fit / norm_t * count / (count - 9));
            if (!(estimate > (trical_wide_t)TRICAL_NOISE_FLOOR)) {
                estimate = (trical_wide_t)TRICAL_NOISE_FLOOR;
            }
            if (!(estimate < (trical_wide_t)TRICAL_NOISE_MAX)) {
                estimate = (trical_wide_t)TRICAL_NOISE_MAX;
            }
        }

        if (estimate > 0 && estimate < _noise_assumed(instance)) {
            instance->priv.noise_estimate = (float)estimate;
        } else if (estimate > (trical_wide_t)TRICAL_NOISE_RAISE *
                              _noise_assumed(instance) &&
                change < (trical_wide_t)1e-8 &&
                overlap > (trical_wide_t)(1.0 - 1e-5)) {
            /*
            The pull of the prior on theta is the gradient of its penalty,
            less the part along theta: g = w (Q theta - (theta' Q theta)
            theta). Were the prior weakened by a factor e, the fit would to
            first order move by (1 - e) C_e^-1 g, where C_e is the matrix
            with the weakened prior, and the quantity minimised would fall
            by (1 - e)^2 g' C_e^-1 g. (The factorisation in `c` is finished
            with, so C_e can go there.)
            */
            const trical_wide_t weaker = (trical_wide_t)TRICAL_NOISE_WEAKER;
            trical_wide_t g[10], y[10], penalty;

            _prior_product(&prior, theta, g);
            penalty = _dot10(theta, g);
            for (i = 0; i < 10; i++) {
                g[i] = weight * (g[i] - penalty * theta[i]);
                y[i] = g[i];
            }

            if (_factorise(instance->priv.moments, xi, &prior, lambda,
                           weaker * weight, c, &shift, &budget) == 0) {
                trical_cholesky_solve(c, 10, y);

                if ((1 - weaker) * (1 - weaker) * _dot10(g, y) <
                        (trical_wide_t)TRICAL_NOISE_MISFIT * fit) {
                    instance->priv.noise_estimate = (float)estimate;
                }
            }
        }
    }

    _gate_update(instance);

    /*
    From now on measurements can be judged by the estimate: the magnitude
    of a calibrated measurement should be within a few times the noise of
    the field norm -- or within a fixed fraction of it, if that is more,
    since the estimate may still be poor where there have been few
    measurements.
    */
    if (_publish(instance, centre, a, 0) == 0 &&
            count >= (trical_wide_t)TRICAL_OUTLIER_MIN_COUNT) {
        const trical_wide_t noise = (trical_wide_t)TRICAL_OUTLIER_NOISE *
                                    _noise_assumed(instance);

        instance->priv.fit_gate = noise > (trical_wide_t)TRICAL_OUTLIER_FIT ?
                                  (float)noise : TRICAL_OUTLIER_FIT;
    }
}

/*
_trical_variance
Writes approximate variances of the elements of `state` to `variance`.

To first order the covariance of the unit quadric vector theta is
    (s^2 theta' N_T theta / count) (M - lambda N)^+
(cf. [1], eq. 67, with the weights replaced by their mean), where ^+ is the
inverse on the complement of theta. This is propagated to the state through
a numerical Jacobian of the quadric-to-calibration conversion.
*/
void _trical_variance(const TRICAL_instance_t *instance,
float variance[TRICAL_STATE_DIM]) {
    const trical_wide_t step = (sizeof(trical_wide_t) == sizeof(float)) ?
        (trical_wide_t)3e-3 : (trical_wide_t)1e-4;
    trical_wide_t c[100];
    trical_wide_t xi[10], u[10], y[10], full[9], base_inverse[9], centre[3];
    trical_wide_t a[9], current[9];
    float jac[TRICAL_STATE_DIM][10];
    const trical_wide_t *theta;
    trical_wide_t count, norm_t, noise, factor, lambda, weight, shift;
    _prior_t prior;
    unsigned int i, j, k, side, budget = TRICAL_FIT_BUDGET;

    assert(instance && variance);

    /*
    Until there is something better, the nominal prior: a tenth of the size
    of the ellipsoid in the bias, and a tenth in the scale factors, in the
    units that the current scale estimate implies (its mean squared
    singular value being the conversion).
    */
    {
        float units2 = 0.0f, bias, scale = 1e-2f;

        for (i = 0; i < 9; i++) {
            const float element = instance->state[3 + i] +
                                  (i % 4 == 0 ? 1.0f : 0.0f);

            units2 += element * element / 3.0f;
        }
        bias = 1e-2f * instance->field_norm * instance->field_norm;
        if (units2 > 0.0f && _float_isfinite(units2) &&
                _float_isfinite(bias / units2)) {
            bias /= units2;
            scale *= units2;
        }
        if (!_float_isfinite(bias)) {
            bias = 1e-2f;
        }

        for (i = 0; i < 3; i++) {
            variance[i] = bias;
        }
        for (i = 3; i < TRICAL_STATE_DIM; i++) {
            variance[i] = scale;
        }
    }

    count = instance->priv.moments[0];
    if (!instance->priv.started || !(count >= 1) ||
            memcmp(instance->state, instance->priv.published,
                   sizeof(instance->state)) != 0) {
        return;
    }

    theta = instance->priv.quadric;
    _trical_basis_mean(instance->priv.moments, xi);
    shift = (trical_wide_t)instance->priv.shift;
    if (_prior_init(instance, xi, &prior) != 0) {
        return;
    }
    _regularisation(instance, xi, &prior, &lambda, &weight);
    if (_factorise(instance->priv.moments, xi, &prior, lambda, weight,
                   c, &shift, &budget) != 0 ||
            _trical_quadric_to_calibration(theta, centre, a,
                                           base_inverse) != 0) {
        return;
    }

    _trical_normalisation_product(xi, theta, 0, y);
    norm_t = _dot10(theta, y);
    noise = _noise_assumed(instance);
    factor = noise * noise * norm_t / count;
    if (!(factor > 0) || !_wide_isfinite(factor)) {
        return;
    }

    /*
    The state is (origin + scale centre, R A - I) with A = gain a and R the
    rotation, held fixed: R A = (I + D) a0^-1 a, where a0 is the current a.
    */
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            current[i * 3 + j] =
                (trical_wide_t)instance->state[3 + i * 3 + j] +
                (i == j ? 1 : 0);
        }
    }
    _mul33(current, base_inverse, full);

    /* numerical Jacobian of the state with respect to theta */
    for (k = 0; k < 10; k++) {
        trical_wide_t value[2][TRICAL_STATE_DIM];

        for (side = 0; side < 2; side++) {
            trical_wide_t perturbed[10], product[9];

            memcpy(perturbed, theta, sizeof(perturbed));
            perturbed[k] += side ? -step : step;

            if (_trical_quadric_to_calibration(perturbed, centre, a,
                                               NULL) != 0) {
                return;
            }
            _mul33(full, a, product);

            for (i = 0; i < 3; i++) {
                value[side][i] = (trical_wide_t)instance->priv.scale *
                                 centre[i];
            }
            for (i = 0; i < 9; i++) {
                value[side][3 + i] = product[i];
            }
        }

        for (i = 0; i < TRICAL_STATE_DIM; i++) {
            const trical_wide_t slope = (value[0][i] - value[1][i]) /
                                        (2 * step);

            if (!_wide_isfloat(slope)) {
                return;
            }
            jac[i][k] = (float)slope;
        }
    }

    /*
    variance[i] = factor * jac[i] V jac[i]', with V = P C^-1 P and
    P = I - theta theta' (which projects theta out)
    */
    for (i = 0; i < TRICAL_STATE_DIM; i++) {
        trical_wide_t along = 0, result;

        for (k = 0; k < 10; k++) {
            along += (trical_wide_t)jac[i][k] * theta[k];
        }
        for (k = 0; k < 10; k++) {
            u[k] = (trical_wide_t)jac[i][k] - along * theta[k];
            y[k] = u[k];
        }
        trical_cholesky_solve(c, 10, y);

        /* P y . u, and u is already orthogonal to theta */
        result = factor * _dot10(u, y);

        /* written so that a NaN leaves the default in place */
        if (result >= 0) {
            variance[i] = _wide_isfloat(result) ? (float)result : FLT_MAX;
        }
    }
}
