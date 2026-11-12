#!/usr/bin/env python
# coding=utf-8

#Copyright (C) 2013 Ben Dyer
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in
#all copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

"""
This is a simple Python interface wrapper around TRICAL, designed more as a
convenient interface for tests than anything else.

If run directly (i.e. `python -m TRICAL`), we expect comma-separated readings
on stdin (3 values per line, ending with \\n), and write calibrated values on
stdout in the same format.

The shared library is looked for, in order: at the path given by the
TRICAL_LIBRARY environment variable; next to this file (where `setup.py`
puts it); and wherever the system looks for shared libraries.
"""

import ctypes
import ctypes.util
import math
import os

__version__ = "2.0.0"

_STATE_DIM = 12

# TRICAL.c checks at compile time that an instance fits in this many bytes,
# so instances can be allocated here without knowing the private part of
# TRICAL_instance_t (which depends on how the library was built).
_INSTANCE_SIZE = 2048

_LIBRARY_NAMES = ("libTRICAL.so", "libTRICAL.dylib", "TRICAL.dll",
                  "libTRICAL.dll")

_TRICAL = None

_Vector = ctypes.c_float * 3
_Matrix = ctypes.c_float * 9


class _Public(ctypes.Structure):
    """
    The public part of TRICAL_instance_t, which is at its start (see
    TRICAL.h).
    """
    _fields_ = [
        ("field_norm", ctypes.c_float),
        ("measurement_noise", ctypes.c_float),
        ("state", ctypes.c_float * _STATE_DIM),
        ("measurement_count", ctypes.c_uint)
    ]


def _find_library():
    """
    Returns the candidate paths (or names) of the TRICAL shared library, most
    preferred first.
    """
    candidates = []

    override = os.environ.get("TRICAL_LIBRARY")
    if override:
        candidates.append(override)

    here = os.path.dirname(os.path.abspath(__file__))
    for name in _LIBRARY_NAMES:
        path = os.path.join(here, name)
        if os.path.exists(path):
            candidates.append(path)

    found = ctypes.util.find_library("TRICAL")
    if found:
        candidates.append(found)

    return candidates


def _init():
    """
    Loads the TRICAL library and sets up the ctypes interface.

    Called automatically the first time a Python instance is created.
    """
    global _TRICAL

    lib = None
    errors = []
    for candidate in _find_library():
        try:
            lib = ctypes.CDLL(candidate)
            break
        except OSError as e:
            errors.append("%s: %s" % (candidate, e))

    if lib is None:
        raise OSError(
            "Could not load the TRICAL shared library. Build it (see "
            "README.md) and either copy it next to %s or set the "
            "TRICAL_LIBRARY environment variable to its path. %s" %
            (os.path.abspath(__file__), "; ".join(errors)))

    if not hasattr(lib, "TRICAL_noise_estimate_get"):
        raise OSError("The TRICAL shared library that was found is older "
                      "than this Python module (%s)" % __version__)

    instance = ctypes.c_void_p

    lib.TRICAL_init.argtypes = [instance]
    lib.TRICAL_init.restype = None

    lib.TRICAL_reset.argtypes = [instance]
    lib.TRICAL_reset.restype = None

    lib.TRICAL_norm_set.argtypes = [instance, ctypes.c_float]
    lib.TRICAL_norm_set.restype = None

    lib.TRICAL_norm_get.argtypes = [instance]
    lib.TRICAL_norm_get.restype = ctypes.c_float

    lib.TRICAL_noise_set.argtypes = [instance, ctypes.c_float]
    lib.TRICAL_noise_set.restype = None

    lib.TRICAL_noise_get.argtypes = [instance]
    lib.TRICAL_noise_get.restype = ctypes.c_float

    lib.TRICAL_noise_estimate_get.argtypes = [instance]
    lib.TRICAL_noise_estimate_get.restype = ctypes.c_float

    lib.TRICAL_measurement_count_get.argtypes = [instance]
    lib.TRICAL_measurement_count_get.restype = ctypes.c_uint

    lib.TRICAL_estimate_update.argtypes = [instance,
                                           ctypes.POINTER(_Vector),
                                           ctypes.POINTER(_Vector)]
    lib.TRICAL_estimate_update.restype = None

    lib.TRICAL_estimate_get.argtypes = [instance,
                                        ctypes.POINTER(_Vector),
                                        ctypes.POINTER(_Matrix)]
    lib.TRICAL_estimate_get.restype = None

    lib.TRICAL_estimate_get_ext.argtypes = [instance,
                                            ctypes.POINTER(_Vector),
                                            ctypes.POINTER(_Matrix),
                                            ctypes.POINTER(_Vector),
                                            ctypes.POINTER(_Matrix)]
    lib.TRICAL_estimate_get_ext.restype = None

    lib.TRICAL_measurement_calibrate.argtypes = [instance,
                                                 ctypes.POINTER(_Vector),
                                                 ctypes.POINTER(_Vector)]
    lib.TRICAL_measurement_calibrate.restype = None

    _TRICAL = lib


def _vector(values, what):
    """
    Converts a sequence of 3 numbers to a C array, checking that it is one.
    """
    try:
        values = tuple(float(v) for v in values)
    except (TypeError, ValueError):
        raise ValueError("%s must be a sequence with 3 items" % what)

    if len(values) != 3:
        raise ValueError("%s must be a sequence with 3 items" % what)

    return _Vector(*values)


class Instance(object):
    def __init__(self, field_norm=1.0, measurement_noise=1e-6):
        """
        Create a new TRICAL instance with the supplied field norm (magnitude)
        and measurement noise.

        The field norm should be the expected magnitude of the field at the
        sensor location, and the measurement noise should be the standard
        deviation of the error in sensor readings (the error being presumed
        to be white Gaussian noise). The latter is only a starting point:
        TRICAL estimates the noise level from the readings themselves (see
        `noise_estimate`).
        """
        if not _TRICAL:
            _init()

        # Sanity-check the input parameters
        if not field_norm > 0.0 or math.isinf(field_norm):
            raise ValueError("Field norm must be > 0.0 (got %f)" % field_norm)

        if not measurement_noise > 0.0 or math.isinf(measurement_noise):
            raise ValueError("Measurement noise must be > 0.0 (got %f)" %
                             measurement_noise)

        # Initialize the internal (C) instance. An array of doubles has the
        # alignment that the instance needs.
        self._buffer = (ctypes.c_double * (_INSTANCE_SIZE // 8))()
        self._pointer = ctypes.addressof(self._buffer)
        self._instance = _Public.from_buffer(self._buffer)

        _TRICAL.TRICAL_init(self._pointer)
        _TRICAL.TRICAL_norm_set(self._pointer, field_norm)
        _TRICAL.TRICAL_noise_set(self._pointer, measurement_noise)

    def __repr__(self):
        return "TRICAL.Instance(%s)" % str({
            "field_norm": self.field_norm,
            "measurement_noise": self.measurement_noise,
            "bias": self.bias,
            "scale": self.scale,
            "measurement_count": self.measurement_count
        })

    def update(self, measurement, reference_field=None):
        """
        Update the calibration estimate based on a new measurement.

        If the expected direction of the field in the frame the calibrated
        readings should be in is known, pass it as `reference_field` (only
        its direction matters); TRICAL then also corrects for the
        misalignment of the sensor. Without it, the scale matrix is
        symmetric.
        """
        measurement = _vector(measurement, "Measurement")
        if reference_field is not None:
            reference_field = _vector(reference_field, "Reference field")

        _TRICAL.TRICAL_estimate_update(self._pointer, measurement,
                                       reference_field)

    def calibrate(self, measurement):
        """
        Given a measurement, return a calibrated measurement based on the
        current calibration estimate.
        """
        calibrated_measurement = _Vector()
        _TRICAL.TRICAL_measurement_calibrate(
            self._pointer, _vector(measurement, "Measurement"),
            calibrated_measurement)

        return tuple(calibrated_measurement)

    def reset(self):
        """
        Discard everything learned so far; the field norm and measurement
        noise are kept.
        """
        _TRICAL.TRICAL_reset(self._pointer)

    @property
    def field_norm(self):
        return _TRICAL.TRICAL_norm_get(self._pointer)

    @field_norm.setter
    def field_norm(self, value):
        if not value > 0.0 or math.isinf(value):
            raise ValueError("Field norm must be > 0.0 (got %f)" % value)
        _TRICAL.TRICAL_norm_set(self._pointer, value)

    @property
    def measurement_noise(self):
        return _TRICAL.TRICAL_noise_get(self._pointer)

    @measurement_noise.setter
    def measurement_noise(self, value):
        if not value > 0.0 or math.isinf(value):
            raise ValueError("Measurement noise must be > 0.0 (got %f)" %
                             value)
        _TRICAL.TRICAL_noise_set(self._pointer, value)

    @property
    def noise_estimate(self):
        """
        The measurement noise (standard deviation) estimated from the
        measurements, or 0.0 if there is no estimate yet.
        """
        return _TRICAL.TRICAL_noise_estimate_get(self._pointer)

    @property
    def measurement_count(self):
        return _TRICAL.TRICAL_measurement_count_get(self._pointer)

    @property
    def bias(self):
        """The bias estimate, a 3-tuple."""
        return tuple(self._instance.state[0:3])

    @property
    def scale(self):
        """
        The scale estimate D, a 9-tuple in row-major order: a calibrated
        measurement is (I + D)(measurement - bias).
        """
        return tuple(self._instance.state[3:12])

    @property
    def variance(self):
        """
        Approximate variances of the elements of the bias and scale
        estimates: a (3-tuple, 9-tuple) pair.
        """
        bias, scale = _Vector(), _Matrix()
        bias_variance, scale_variance = _Vector(), _Matrix()

        _TRICAL.TRICAL_estimate_get_ext(self._pointer, bias, scale,
                                        bias_variance, scale_variance)
        return tuple(bias_variance), tuple(scale_variance)

    @property
    def state(self):
        """
        The calibration estimate as 12 values: the bias followed by the scale
        matrix. Assigning to this (e.g. a calibration saved earlier) makes
        TRICAL carry on from the values assigned.
        """
        return tuple(self._instance.state)

    @state.setter
    def state(self, values):
        values = tuple(float(v) for v in values)
        if len(values) != _STATE_DIM:
            raise ValueError("State must be a sequence with %d items" %
                             _STATE_DIM)
        self._instance.state[:] = values


def _squared_norm(v):
    return v[0] * v[0] + v[1] * v[1] + v[2] * v[2]


def generate_html_viz(instance, samples):
    """
    Generate a WebGL visualisation of the raw samples, and of the same
    samples calibrated with the estimate as it stood when each arrived.

    Each sample is either a measurement, or a (measurement, reference field)
    pair.
    """
    raw_data = []
    iteratively_calibrated_data = []
    squared_magnitude = 0.0

    # Iterate over the samples and generate the iteratively calibrated data
    for sample in samples:
        if len(sample) == 2:
            measurement, reference_field = sample
        else:
            measurement, reference_field = sample, None

        instance.update(measurement, reference_field)
        calibrated_measurement = instance.calibrate(measurement)
        raw_data.append(measurement)
        iteratively_calibrated_data.append(calibrated_measurement)

        # Update the maximum magnitude seen
        squared_magnitude = max(squared_magnitude, _squared_norm(measurement))
        squared_magnitude = max(squared_magnitude,
                                _squared_norm(calibrated_measurement))

    # And now insert them into the appropriate part of the HTML
    html_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "viz.html")
    with open(html_path, "r") as f:
        html = f.read()

    # Insert the visualisation data -- output measurements as a single array
    # of points
    raw_points = ",\n".join(",".join("%.4f" % v for v in measurement)
                            for measurement in raw_data)
    calibrated_points = ",\n".join(
        ",".join("%.4f" % v for v in measurement)
        for measurement in iteratively_calibrated_data)

    html = html.replace("{{raw}}", "[" + raw_points + "]")
    html = html.replace("{{calibrated}}", "[" + calibrated_points + "]")

    # TODO: maybe insert the calibration estimate as well?
    html = html.replace("{{magnitude}}",
                        "%.3f" % math.sqrt(squared_magnitude))
    html = html.replace("{{fieldNorm}}", "%.3f" % instance.field_norm)

    return html
