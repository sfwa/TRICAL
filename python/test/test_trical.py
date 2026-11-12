"""
Tests of the Python interface. Run with the TRICAL_LIBRARY environment
variable set to the path of the shared library, unless it has been installed
alongside the module:

    TRICAL_LIBRARY=/path/to/libTRICAL.so \\
        python -m unittest discover python/test
"""

import math
import os
import random
import subprocess
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                ".."))

import TRICAL


def _direction(rng):
    while True:
        v = [rng.gauss(0.0, 1.0) for _ in range(3)]
        n = math.sqrt(sum(x * x for x in v))
        if n > 1e-3:
            return [x / n for x in v]


def _readings(rng, count, norm=1.0, noise=0.005):
    """
    Readings of a sensor with a bias and scale errors: (reading, field)
    pairs, where reading = field / scale + bias + noise.
    """
    bias = (0.3 * norm, -0.2 * norm, 0.1 * norm)
    scale = (1.1, 0.9, 1.05)
    out = []
    for _ in range(count):
        field = [norm * x for x in _direction(rng)]
        reading = tuple(field[i] / scale[i] + bias[i] +
                        rng.gauss(0.0, noise * norm) for i in range(3))
        out.append((reading, tuple(field)))
    return out


def _magnitude_error(instance, readings, norm=1.0):
    total = 0.0
    for reading, _ in readings:
        c = instance.calibrate(reading)
        total += (math.sqrt(sum(x * x for x in c)) - norm) ** 2
    return math.sqrt(total / len(readings)) / norm


class TestInstance(unittest.TestCase):
    def test_defaults(self):
        instance = TRICAL.Instance()
        self.assertAlmostEqual(1.0, instance.field_norm)
        self.assertAlmostEqual(1e-6, instance.measurement_noise)
        self.assertEqual((0.0, 0.0, 0.0), instance.bias)
        self.assertEqual((0.0,) * 9, instance.scale)
        self.assertEqual(0, instance.measurement_count)
        self.assertEqual(0.0, instance.noise_estimate)

    def test_invalid_arguments(self):
        self.assertRaises(ValueError, TRICAL.Instance, field_norm=0.0)
        self.assertRaises(ValueError, TRICAL.Instance, measurement_noise=-1.0)

        instance = TRICAL.Instance()
        self.assertRaises(ValueError, instance.update, (1.0, 2.0))
        self.assertRaises(ValueError, instance.update, None)
        self.assertRaises(ValueError, instance.update, (1.0, 2.0, 3.0),
                          (1.0, 2.0))
        self.assertRaises(ValueError, instance.calibrate, "abc")

    def test_uncalibrated_is_identity(self):
        instance = TRICAL.Instance()
        self.assertEqual((1.0, 2.0, 3.0), instance.calibrate((1.0, 2.0, 3.0)))

    def test_calibration(self):
        rng = random.Random(1)
        readings = _readings(rng, 1500)
        instance = TRICAL.Instance(field_norm=1.0, measurement_noise=0.005)

        self.assertGreater(_magnitude_error(instance, readings), 0.1)
        for reading, _ in readings:
            instance.update(reading)

        self.assertEqual(1500, instance.measurement_count)
        self.assertLess(_magnitude_error(instance, readings), 0.01)
        self.assertAlmostEqual(0.3, instance.bias[0], delta=0.01)
        self.assertAlmostEqual(-0.2, instance.bias[1], delta=0.01)
        self.assertAlmostEqual(0.1, instance.bias[2], delta=0.01)
        self.assertAlmostEqual(0.1, instance.scale[0], delta=0.01)
        self.assertAlmostEqual(-0.1, instance.scale[4], delta=0.01)
        self.assertAlmostEqual(0.05, instance.scale[8], delta=0.01)
        self.assertAlmostEqual(0.005, instance.noise_estimate, delta=0.001)

        bias_variance, scale_variance = instance.variance
        self.assertEqual(3, len(bias_variance))
        self.assertEqual(9, len(scale_variance))
        self.assertTrue(all(0.0 <= v < 1e-4 for v in bias_variance))

    def test_reference_field(self):
        rng = random.Random(2)
        readings = _readings(rng, 1500, norm=50.0)
        instance = TRICAL.Instance(field_norm=50.0, measurement_noise=0.25)

        for reading, field in readings:
            instance.update(reading, field)

        worst = 0.0
        for reading, field in readings[:100]:
            c = instance.calibrate(reading)
            worst = max(worst, max(abs(c[i] - field[i]) for i in range(3)))
        self.assertLess(worst, 50.0 * 0.03)

    def test_reset_and_state(self):
        rng = random.Random(3)
        readings = _readings(rng, 1000)
        instance = TRICAL.Instance()
        for reading, _ in readings:
            instance.update(reading)

        saved = instance.state
        self.assertEqual(12, len(saved))

        instance.reset()
        self.assertEqual((0.0,) * 12, instance.state)
        self.assertGreater(_magnitude_error(instance, readings), 0.1)

        # restoring a saved calibration takes effect immediately
        instance.state = saved
        self.assertLess(_magnitude_error(instance, readings), 0.01)
        instance.update(readings[0][0])
        self.assertLess(_magnitude_error(instance, readings), 0.01)

        self.assertRaises(ValueError, setattr, instance, "state", (1.0, 2.0))

    def test_instances_are_independent(self):
        rng = random.Random(4)
        a, b = TRICAL.Instance(), TRICAL.Instance()
        for reading, _ in _readings(rng, 200):
            a.update(reading)
        self.assertEqual((0.0,) * 12, b.state)
        self.assertNotEqual((0.0,) * 12, a.state)

    def test_html(self):
        rng = random.Random(5)
        samples = [reading for reading, _ in _readings(rng, 50)]
        html = TRICAL.generate_html_viz(TRICAL.Instance(), samples)
        self.assertIn("TRICAL.data.raw = [", html)
        self.assertNotIn("{{", html)


class TestCommandLine(unittest.TestCase):
    def _run(self, arguments, text):
        environment = dict(os.environ)
        environment["PYTHONPATH"] = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "..")
        process = subprocess.Popen(
            [sys.executable, "-m", "TRICAL"] + arguments,
            stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, env=environment, universal_newlines=True)
        out, err = process.communicate(text)
        return process.returncode, out, err

    def test_usage(self):
        code, _, err = self._run([], "")
        self.assertEqual(1, code)
        self.assertIn("Usage", err)

    def test_calibration(self):
        rng = random.Random(6)
        readings = _readings(rng, 800)
        text = "".join("%.6f,%.6f,%.6f\n" % r for r, _ in readings)
        text += "not,a,number\n\n1.0,2.0\n"

        code, out, err = self._run(["1.0", "0.005"], text)
        self.assertEqual(0, code)
        self.assertIn("CALIBRATION", err)

        lines = out.strip().split("\n")
        self.assertEqual(800, len(lines))
        last = [float(v) for v in lines[-1].split(",")]
        self.assertAlmostEqual(1.0, math.sqrt(sum(x * x for x in last)),
                               delta=0.03)

    def test_reference_columns(self):
        rng = random.Random(7)
        readings = _readings(rng, 800)
        text = "".join("%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n" % (r + f)
                       for r, f in readings)

        code, out, _ = self._run(["1.0", "0.005"], text)
        self.assertEqual(0, code)
        self.assertEqual(800, len(out.strip().split("\n")))


if __name__ == "__main__":
    unittest.main()
