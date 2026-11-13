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

#ifndef TRICAL_FIT_H
#define TRICAL_FIT_H

#include "TRICAL.h"
#include "linalg.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
Tuning constants. None of these needs to be adjusted for a particular sensor;
they are collected here so that their roles are documented in one place.
Each can be overridden when the library is built (e.g. -DTRICAL_GATE=0.05f).
*/

/*
Measurements are grouped into cells: a new cell starts when a measurement is
further than TRICAL_GATE from the centre of the current cell, as a fraction
of the size of the initial estimate's ellipsoid (i.e. of the field norm in
the units of the measurements) -- or further than TRICAL_GATE_NOISE times
the measurement noise, if that is more, so that noise alone does not start
new cells. The centre of a cell is the mean of its first
TRICAL_GATE_AVERAGE measurements.

The fit is refreshed at the start of each cell, and at most
TRICAL_CELL_QUOTA measurements per cell are used (in proportion more for the
larger cells of a noisy sensor), so that a sensor at rest neither costs
anything nor comes to dominate the fit. (At a typical sample rate a slowly
turning sensor produces a few dozen measurements per cell, all of which are
worth having.)
*/
#ifndef TRICAL_GATE
#define TRICAL_GATE (0.1f)
#endif
#ifndef TRICAL_GATE_NOISE
#define TRICAL_GATE_NOISE (8.0f)
#endif
#ifndef TRICAL_GATE_AVERAGE
#define TRICAL_GATE_AVERAGE (8u)
#endif
#ifndef TRICAL_CELL_QUOTA
#define TRICAL_CELL_QUOTA (64u)
#endif

/*
The fit is refreshed when a measurement starts a new cell, if at least
TRICAL_REFIT_INTERVAL measurements have been used since the last time, or
at least 1 / TRICAL_REFIT_FRACTION of the total if that is fewer (so that
the first few dozen measurements each count at once). That keeps the cost
down when the sensor turns so fast that every measurement starts a cell.

There is no fit at all until TRICAL_REFIT_MIN_COUNT measurements have been
used: fewer than there are parameters can only add their noise to the
initial estimate.
*/
#ifndef TRICAL_REFIT_INTERVAL
#define TRICAL_REFIT_INTERVAL (8u)
#endif
#ifndef TRICAL_REFIT_FRACTION
#define TRICAL_REFIT_FRACTION (16u)
#endif
#ifndef TRICAL_REFIT_MIN_COUNT
#define TRICAL_REFIT_MIN_COUNT (9u)
#endif

/*
With single-precision sums, the sums are halved when this many measurements
have been used, so that they stay within the precision available (see
_trical_refit). Zero means never, which is right for double precision.
*/
#ifndef TRICAL_COUNT_MAX
#define TRICAL_COUNT_MAX ((sizeof(trical_wide_t) == sizeof(float)) ? \
                          131072.0f : 0.0f)
#endif

/*
Strength of the prior that holds the estimate at its initial value where the
measurements do not determine it. With W measurements used and measurement
noise of standard deviation s (as a fraction of the field norm), the prior
weighs as much as s^2 (TRICAL_PRIOR_BASE + TRICAL_PRIOR_GROWTH sqrt(W))
measurements.

The first term corresponds to an initial uncertainty of about a tenth of the
field norm in the bias and a tenth in the scale factors. The second keeps the
prior above the sampling fluctuations of the noise correction (which grow as
sqrt(W)), so that directions that are never observed stay put indefinitely,
while its weight relative to the data still vanishes as 1 / sqrt(W).

The scale factors are held to their initial values in that way however far
off those turn out to be. The bias is not: the hold on it weakens once the
measurements are, in the RMS sense, further than about
TRICAL_PRIOR_BIAS_RANGE field norms from the ellipsoid of the initial
estimate, so that an initial bias that is wrong by several times the field
norm does not stand in the way.
*/
#ifndef TRICAL_PRIOR_BASE
#define TRICAL_PRIOR_BASE (300.0f)
#endif
#ifndef TRICAL_PRIOR_GROWTH
#define TRICAL_PRIOR_GROWTH (20.0f)
#endif
#ifndef TRICAL_PRIOR_BIAS_RANGE
#define TRICAL_PRIOR_BIAS_RANGE (0.5f)
#endif

/*
That is for an initial estimate that the measurements show to be wrong, which
must not get in their way. One that they do not contradict can be held much
more firmly, and should be: TRICAL_PRIOR_HOLD is added to
TRICAL_PRIOR_GROWTH while the ellipsoid of the initial estimate fits the
measurements as well as the best fit does, and progressively less of it as
the mean squared distance of the measurements from that ellipsoid exceeds
that of the best fit, by a fraction of the order of
TRICAL_PRIOR_HOLD_EXCESS. That is what keeps a good calibration from
drifting while the sensor is only turning about one axis, say.
*/
#ifndef TRICAL_PRIOR_HOLD
#define TRICAL_PRIOR_HOLD (400.0f)
#endif
#ifndef TRICAL_PRIOR_HOLD_EXCESS
#define TRICAL_PRIOR_HOLD_EXCESS (0.25f)
#endif

/*
The running sums are only as good as the arithmetic: after W additions the
rounding error of a sum, relative to its size, is of the order of
sqrt(W) times the precision of the type. Whatever the measurement noise, the
prior is never weaker than that, or rounding errors would pass for evidence
in the directions that the measurements say little about; and never weaker
than TRICAL_PRIOR_ROUNDING times that while the measurements do not
contradict the initial estimate (as for TRICAL_PRIOR_HOLD). (With
double-precision sums this never comes into play.)
*/
#ifndef TRICAL_PRIOR_ROUNDING
#define TRICAL_PRIOR_ROUNDING (8.0f)
#endif

/*
Likewise for the test of whether the measurements contradict the initial
estimate: a mean squared distance from an ellipsoid is the small difference
of large sums, and cannot be told from zero if it is less than
TRICAL_PRIOR_HOLD_ROUNDING times the rounding error above.
*/
#ifndef TRICAL_PRIOR_HOLD_ROUNDING
#define TRICAL_PRIOR_HOLD_ROUNDING (64.0f)
#endif

/*
The hold on the overall scale factor is formulated about the centre of the
previous fit. If the centre of a fit is further than TRICAL_PRIOR_RECENTRE
field norms from where it was taken to be, the fit is repeated about the new
centre, at most TRICAL_PRIOR_RECENTRE_LIMIT times in all.
*/
#ifndef TRICAL_PRIOR_RECENTRE
#define TRICAL_PRIOR_RECENTRE (0.2f)
#endif
#ifndef TRICAL_PRIOR_RECENTRE_LIMIT
#define TRICAL_PRIOR_RECENTRE_LIMIT (4u)
#endif

/*
Until the fit yields an estimate of the measurement noise, the value set by
the user is assumed, limited to this range (as fractions of the field norm).
The lower limit keeps the prior effective for the first few measurements
whatever the user's setting; the upper limit keeps it from overwhelming them.
*/
#ifndef TRICAL_NOISE_INITIAL_MIN
#define TRICAL_NOISE_INITIAL_MIN (0.03f)
#endif
#ifndef TRICAL_NOISE_INITIAL_MAX
#define TRICAL_NOISE_INITIAL_MAX (0.05f)
#endif

/*
Minimum number of measurements before the noise is estimated, and the
smallest estimate (as a fraction of the field norm): below this the
arithmetic of the running sums, not the sensor, is what is being measured.
*/
#ifndef TRICAL_NOISE_MIN_COUNT
#define TRICAL_NOISE_MIN_COUNT (16)
#endif
#ifndef TRICAL_NOISE_FLOOR
#define TRICAL_NOISE_FLOOR (1e-4f)
#endif

/*
Upper limit on the noise variance that the fit corrects for, as a fraction
of the squared size of the ellipsoid. And if a fit turns out to be closer to
the measurements than TRICAL_CORRECTION_REPEAT times the noise variance it
was corrected for, that was too much (the previous fit was a poor guide to
this one), and the fit is repeated with a correction for what it found.
*/
#ifndef TRICAL_CORRECTION_MAX
#define TRICAL_CORRECTION_MAX (1e-2f)
#endif
#ifndef TRICAL_CORRECTION_REPEAT
#define TRICAL_CORRECTION_REPEAT (0.5f)
#endif

/*
The noise estimate is raised when the fit suggests more than
TRICAL_NOISE_RAISE times the current one (less is not worth the test that
follows), but only when the fit is determined by the data rather than held
back by the prior: when weakening the prior by a factor of
TRICAL_NOISE_WEAKER would improve the fit by no more than the fraction
TRICAL_NOISE_MISFIT of the squared distance of the measurements from the
ellipsoid. It is never raised beyond TRICAL_NOISE_MAX (as a fraction of the
field norm): a calibration cannot be expected of measurements noisier than
that.
*/
#ifndef TRICAL_NOISE_RAISE
#define TRICAL_NOISE_RAISE (1.05f)
#endif
#ifndef TRICAL_NOISE_WEAKER
#define TRICAL_NOISE_WEAKER (0.1f)
#endif
#ifndef TRICAL_NOISE_MISFIT
#define TRICAL_NOISE_MISFIT (0.1f)
#endif
#ifndef TRICAL_NOISE_MAX
#define TRICAL_NOISE_MAX (0.1f)
#endif

/*
Weight of the prior on the rotation, relative to the reference data: only
there to define the rotation about axes that the data leave undetermined.
*/
#ifndef TRICAL_ROTATION_PRIOR
#define TRICAL_ROTATION_PRIOR (1e-3f)
#endif

/*
The origin of the normalised coordinates is moved to the mean of the data
when that is further than this from the current origin (in field norms),
which keeps the fit well conditioned.
*/
#ifndef TRICAL_RECENTRE
#define TRICAL_RECENTRE (1.5f)
#endif

/*
Measurements further than this (in field norms) from the initial bias
estimate are ignored until there are enough measurements to judge by; after
that, measurements further than TRICAL_OUTLIER_RANGE times (1 + the RMS
radius of the data) from the mean of the data are ignored.
*/
#ifndef TRICAL_RANGE
#define TRICAL_RANGE (8.0f)
#endif
#ifndef TRICAL_OUTLIER_RANGE
#define TRICAL_OUTLIER_RANGE (4.0f)
#endif
#ifndef TRICAL_OUTLIER_MIN_COUNT
#define TRICAL_OUTLIER_MIN_COUNT (16)
#endif

/*
Once there is an estimate based on at least TRICAL_OUTLIER_MIN_COUNT
measurements, a measurement is also ignored if the magnitude of its
calibrated value differs from the field norm by more than
TRICAL_OUTLIER_NOISE times the measurement noise, or by the fraction
TRICAL_OUTLIER_FIT of the field norm if that is more. (The estimate may
well be that far out where there have been no measurements yet. A single
stray measurement, on the other hand, has an influence that grows as the
fourth power of its distance, so it does not have to be very far out to do
damage.)

That holds for at most TRICAL_OUTLIER_RUN such measurements in a row. If
there are more, and they do not all come from one place, they are used after
all: the estimate is then what is wrong.
*/
#ifndef TRICAL_OUTLIER_NOISE
#define TRICAL_OUTLIER_NOISE (8.0f)
#endif
#ifndef TRICAL_OUTLIER_FIT
#define TRICAL_OUTLIER_FIT (0.25f)
#endif
#ifndef TRICAL_OUTLIER_RUN
#define TRICAL_OUTLIER_RUN (64u)
#endif

/*
If this many consecutive measurements are ignored as out of range, and they
do not all come from one place (as those of a stuck or saturated sensor
do), the sensor has evidently changed, and estimation restarts where the
measurements now are, keeping only the scale factors of the current
estimate.
*/
#ifndef TRICAL_REJECT_LIMIT
#define TRICAL_REJECT_LIMIT (256u)
#endif

/*
The range of sizes that the ellipsoid of the measurements may have, in the
units of the measurements: within it, the squares of distances between
measurements are safely representable in single precision.
*/
#ifndef TRICAL_SCALE_MIN
#define TRICAL_SCALE_MIN (1e-15f)
#endif
#ifndef TRICAL_SCALE_MAX
#define TRICAL_SCALE_MAX (1e15f)
#endif

/*
Limits on what is accepted as a valid fit: the largest ratio between the
principal scale factors; and the smallest and largest semi-axis, relative to
the size of the ellipsoid of the initial estimate. The largest is the
distance from the initial bias at which measurements are still accepted
(TRICAL_RANGE): a sensor whose measurements are used at all can then be
learned. The smallest is about the size of a cell of measurements
(TRICAL_GATE), below which a moving sensor cannot be told from one at rest
anyway.
*/
#ifndef TRICAL_AXIS_RATIO_MAX
#define TRICAL_AXIS_RATIO_MAX (8.0f)
#endif
#ifndef TRICAL_SIZE_MIN
#define TRICAL_SIZE_MIN (0.125f)
#endif
#ifndef TRICAL_SIZE_MAX
#define TRICAL_SIZE_MAX (8.0f)
#endif

/*
If a fit is not an acceptable ellipsoid it is repeated, at most
TRICAL_FIT_ATTEMPTS times: first without the noise correction, then with the
hold on the scale factors strengthened by a factor of TRICAL_FIT_PRIOR_STEP
each time.
*/
#ifndef TRICAL_FIT_ATTEMPTS
#define TRICAL_FIT_ATTEMPTS (5u)
#endif
#ifndef TRICAL_FIT_PRIOR_STEP
#define TRICAL_FIT_PRIOR_STEP (8.0f)
#endif

/*
The most matrix factorisations (the bulk of the work) that may be made for
one refit, all told. A refit normally takes one; the repeats and fallbacks
above take more, but only for measurements that determine very little or
are far from where the previous fit put them. If the limit is reached the
estimate is left as it was.
*/
#ifndef TRICAL_FIT_BUDGET
#define TRICAL_FIT_BUDGET (24u)
#endif

/*
Weight of the hold on the shape (in the units of the matrix of the
eigenproblem) with which the best fit of an ellipsoid of the initial shape
is found: enough to leave the measurements no say in the shape.
*/
#ifndef TRICAL_FIT_PROBE
#define TRICAL_FIT_PROBE (1e4f)
#endif

/*
Inverse iteration. The matrix is shifted to make it safely positive definite
for the Cholesky factorisation: by at least TRICAL_SHIFT (relative to
normalised coordinates, in which its significant eigenvalues are of the
order of 0.01 to 10), and by ten times more at each further attempt, up to
TRICAL_SHIFT_MAX. Each fit starts a step below the shift that the last one
needed.
*/
#ifndef TRICAL_SHIFT
#define TRICAL_SHIFT ((sizeof(trical_wide_t) == sizeof(float)) ? \
                      (trical_wide_t)1e-4 : (trical_wide_t)1e-6)
#endif
#ifndef TRICAL_SHIFT_MAX
#define TRICAL_SHIFT_MAX (1.0f)
#endif
#ifndef TRICAL_ITERATIONS_MIN
#define TRICAL_ITERATIONS_MIN (3u)
#endif
#ifndef TRICAL_ITERATIONS_MAX
#define TRICAL_ITERATIONS_MAX (12u)
#endif

/* Internal function prototypes */

/*
_trical_moments_accumulate
Adds the 35 monomials of degree <= 4 of the normalised measurement `x` to
`sums`. The products are formed in single precision.
*/
void _trical_moments_accumulate(trical_wide_t sums[TRICAL_NUM_MOMENTS],
const float x[3]);

/*
_trical_moments_shift
Converts sums of monomials of x to sums of monomials of (x - delta).
*/
void _trical_moments_shift(trical_wide_t sums[TRICAL_NUM_MOMENTS],
const trical_wide_t delta[3]);

/*
_trical_basis_mean
Writes to `xi` the mean of the quadric basis
xi = (xx, yy, zz, 2yz, 2xz, 2xy, 2x, 2y, 2z, 1) over the measurements, from
the monomial sums.
*/
void _trical_basis_mean(const trical_wide_t sums[TRICAL_NUM_MOMENTS],
trical_wide_t xi[10]);

/*
_trical_moment_product
y = M theta, where M = E[xi xi^T] is the mean moment matrix of the quadric
basis.
*/
void _trical_moment_product(const trical_wide_t sums[TRICAL_NUM_MOMENTS],
const trical_wide_t theta[10], trical_wide_t y[10]);

/*
_trical_normalisation_product
y = N theta, where N is the Taubin normalisation matrix N_T = E[J J^T] (J
being the Jacobian of xi with respect to x) for the measurements whose mean
quadric basis is `xi` -- plus, if `hyper` is non-zero, the second-order
("hyper") correction 2 S[xi e^T].
*/
void _trical_normalisation_product(const trical_wide_t xi[10],
const trical_wide_t theta[10], int hyper, trical_wide_t y[10]);

/*
_trical_quadric_from_calibration
Writes to `quadric` the unit-norm coefficients of the ellipsoid
|A (x - c)| = 1, where `ata` is A^T A.
*/
void _trical_quadric_from_calibration(const trical_wide_t ata[9],
const trical_wide_t centre[3], trical_wide_t quadric[10]);

/*
_trical_quadric_to_calibration
The inverse of the above, with A symmetric positive definite. Returns 0 on
success, or -1 if `quadric` is not an ellipsoid of acceptable shape.
`a_inverse` may be NULL.
*/
int _trical_quadric_to_calibration(const trical_wide_t quadric[10],
trical_wide_t centre[3], trical_wide_t a[9], trical_wide_t a_inverse[9]);

/*
_trical_quadric_translate
Re-expresses `quadric` in coordinates x' = x - delta, and normalises it.
*/
void _trical_quadric_translate(trical_wide_t quadric[10],
const trical_wide_t delta[3]);

/*
_trical_polar_rotation
Writes to `q` the rotation matrix closest to `k` (its orthogonal polar
factor, with the determinant forced to +1). Returns 0 on success, or -1 if
`k` is too close to zero.
*/
int _trical_polar_rotation(const trical_wide_t k[9], trical_wide_t q[9]);

/*
_trical_recentre
Moves the origin of the normalised coordinates by (approximately) `shift`,
re-expressing everything that is stored in those coordinates.
*/
void _trical_recentre(TRICAL_instance_t *instance,
const trical_wide_t shift[3]);

/*
_trical_start
Starts accumulating, taking the current `state` of `instance` as the initial
estimate. Returns 0 on success, or -1 if `state` is not a usable calibration.
*/
int _trical_start(TRICAL_instance_t *instance, const float bias[3]);

/*
_trical_republish
Writes to `state` the calibration that the current solution of `instance`
amounts to with its current field norm. Returns 0 on success, or -1 (and
`state` is untouched) if that is not a set of finite numbers.
*/
int _trical_republish(TRICAL_instance_t *instance);

/*
_trical_refit
Refits the calibration to the measurements accumulated so far, and updates
`state` if a valid fit is obtained.
*/
void _trical_refit(TRICAL_instance_t *instance);

/*
_trical_variance
Writes approximate variances of the elements of `state` to `variance`.
*/
void _trical_variance(const TRICAL_instance_t *instance,
float variance[TRICAL_STATE_DIM]);

#ifdef __cplusplus
}
#endif

#endif
