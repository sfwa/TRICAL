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
Timing of the public API on a simulated sensor: mean and worst-case time per
call of TRICAL_estimate_update for a tumbling sensor (at several rates of
rotation per sample) and for a sensor at rest, and of
TRICAL_measurement_calibrate.

Only the public API is used, so this also builds against other versions of
the library.
*/

#define _POSIX_C_SOURCE 199309L

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "TRICAL.h"

#define SAMPLES 200000

static float readings[SAMPLES][3];
static float references[SAMPLES][3];

static double now(void) {
    struct timespec t;

    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec * 1e9 + (double)t.tv_nsec;
}

/* xorshift: good enough for noise */
static unsigned int rng_state = 2463534242u;

static double uniform(void) {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 17;
    rng_state ^= rng_state << 5;
    return ((double)rng_state + 0.5) / 4294967296.0;
}

static double normal(void) {
    return sqrt(-2.0 * log(uniform())) * cos(6.283185307179586 * uniform());
}

/*
Readings of a sensor with a bias and scale errors, in a field of magnitude
1, turning by `step` radians per sample about a slowly wandering axis.
*/
static void simulate(double step, double noise) {
    double f[3] = { 1.0, 0.0, 0.0 }, axis[3] = { 0.0, 0.0, 1.0 };
    int i, k;

    for (i = 0; i < SAMPLES; i++) {
        double n, c = cos(step), s = sin(step), dot, cross[3], h[3];

        for (k = 0; k < 3; k++) {
            axis[k] += 0.05 * normal();
        }
        n = sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
        for (k = 0; k < 3; k++) {
            axis[k] /= n;
        }

        /* Rodrigues' rotation of f about axis */
        dot = axis[0] * f[0] + axis[1] * f[1] + axis[2] * f[2];
        cross[0] = axis[1] * f[2] - axis[2] * f[1];
        cross[1] = axis[2] * f[0] - axis[0] * f[2];
        cross[2] = axis[0] * f[1] - axis[1] * f[0];
        for (k = 0; k < 3; k++) {
            h[k] = f[k] * c + cross[k] * s + axis[k] * dot * (1.0 - c);
        }
        n = sqrt(h[0] * h[0] + h[1] * h[1] + h[2] * h[2]);
        for (k = 0; k < 3; k++) {
            f[k] = h[k] / n;
            references[i][k] = (float)f[k];
        }

        /* scale errors of a few percent, bias of 0.3 */
        readings[i][0] = (float)(1.05 * f[0] + 0.02 * f[1] + 0.20 +
                                 noise * normal());
        readings[i][1] = (float)(0.02 * f[0] + 0.97 * f[1] - 0.15 +
                                 noise * normal());
        readings[i][2] = (float)(1.02 * f[2] + 0.16 + noise * normal());
    }
}

static void run(const char *label, int with_reference) {
    TRICAL_instance_t cal;
    double start, total, worst = 0.0, previous;
    float out[3], sink = 0.0f;
    int i;

    TRICAL_init(&cal);
    TRICAL_noise_set(&cal, 0.01f);

    start = now();
    previous = start;
    for (i = 0; i < SAMPLES; i++) {
        double t;

        TRICAL_estimate_update(&cal, readings[i],
                               with_reference ? references[i] : readings[i]);
        t = now();
        if (t - previous > worst) {
            worst = t - previous;
        }
        previous = t;
    }
    total = previous - start;

    TRICAL_measurement_calibrate(&cal, readings[SAMPLES - 1], out);
    sink = out[0] * out[0] + out[1] * out[1] + out[2] * out[2];

    printf("%-34s %8.1f ns/update   worst %8.0f ns   |calibrated| = %.4f\n",
           label, total / SAMPLES, worst, sqrt((double)sink));
}

int main(void) {
    static const double degree = 0.017453292519943295;
    TRICAL_instance_t cal;
    double start, total;
    float out[3], sink = 0.0f;
    int i, repeat;

    printf("sizeof(TRICAL_instance_t) = %u bytes\n",
           (unsigned int)sizeof(TRICAL_instance_t));
    printf("(timings include about 20-40 ns of clock overhead per call)\n");

    simulate(10.0 * degree, 0.01);
    run("tumbling, 10 deg/sample", 0);
    run("tumbling, 10 deg/sample, reference", 1);

    simulate(1.0 * degree, 0.01);
    run("tumbling, 1 deg/sample", 0);
    run("tumbling, 1 deg/sample, reference", 1);

    simulate(0.1 * degree, 0.01);
    run("tumbling, 0.1 deg/sample", 0);
    run("tumbling, 0.1 deg/sample, reference", 1);

    simulate(0.0, 0.01);
    run("at rest", 0);
    run("at rest, reference", 1);

    /* calibration of a reading, without the clock in the loop */
    simulate(1.0 * degree, 0.01);
    TRICAL_init(&cal);
    for (i = 0; i < 5000; i++) {
        TRICAL_estimate_update(&cal, readings[i], readings[i]);
    }
    start = now();
    for (repeat = 0; repeat < 50; repeat++) {
        for (i = 0; i < SAMPLES; i++) {
            TRICAL_measurement_calibrate(&cal, readings[i], out);
            sink += out[0];
        }
    }
    total = now() - start;
    printf("%-34s %8.2f ns/call     (checksum %g)\n",
           "TRICAL_measurement_calibrate", total / (50.0 * SAMPLES),
           (double)sink);

    return 0;
}
