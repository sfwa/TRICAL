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

import os
import sys
from ctypes import *


"""
This is a simple Python interface wrapper around TRICAL, designed more as a
convenient interface for tests than anything else.

If run directly (i.e. `python -m TRICAL`), we read expect comma-separated
readings on stdin (3 values per line, ending with \\n), and write calibrated
values on stdout in the same format.
"""


_TRICAL = None


class _Instance(Structure):
    def __repr__(self):
        fields = {
            "field_norm": self.field_norm,
            "measurement_noise": self.measurement_noise,
            "state": tuple(self.state),
            "state_covariance": tuple(self.state_covariance),
            "measurement_count": self.measurement_count
        }
        return str(fields)


def _init():
    """
    Loads the TRICAL library and sets up the ctypes interface.

    Called automatically the first time a Python instance is created.
    """

    global _TRICAL

    # TODO: search properly
    lib = os.path.join(os.path.dirname(__file__), "libTRICAL.dylib")
    _TRICAL = cdll.LoadLibrary(lib)

    # Set up the _Instance structure based on the definition in TRICAL.h
    _Instance._fields_ = [
        ("field_norm", c_float),
        ("measurement_noise", c_float),
        ("state", c_float * 9),
        ("state_covariance", c_float * 9 * 9),
        ("measurement_count", c_uint)
    ]

    _TRICAL.TRICAL_init.argtypes = [POINTER(_Instance)]
    _TRICAL.TRICAL_init.restype = None

    _TRICAL.TRICAL_norm_set.argtypes = [POINTER(_Instance), c_float]
    _TRICAL.TRICAL_norm_set.restype = None

    _TRICAL.TRICAL_norm_get.argtypes = [POINTER(_Instance)]
    _TRICAL.TRICAL_norm_get.restype = c_float

    _TRICAL.TRICAL_noise_set.argtypes = [POINTER(_Instance), c_float]
    _TRICAL.TRICAL_noise_set.restype = None

    _TRICAL.TRICAL_noise_get.argtypes = [POINTER(_Instance)]
    _TRICAL.TRICAL_noise_get.restype = c_float

    _TRICAL.TRICAL_measurement_count_get.argtypes = [POINTER(_Instance)]
    _TRICAL.TRICAL_measurement_count_get.restype = c_uint

    _TRICAL.TRICAL_estimate_update.argtypes = [POINTER(_Instance),
                                               POINTER(c_float * 3)]
    _TRICAL.TRICAL_estimate_update.restype = None

    _TRICAL.TRICAL_estimate_get.argtypes = [POINTER(_Instance),
                                            POINTER(c_float * 3),
                                            POINTER(c_float * 9)]
    _TRICAL.TRICAL_estimate_get.restype = None

    _TRICAL.TRICAL_estimate_get_ext.argtypes = [POINTER(_Instance),
                                                POINTER(c_float * 3),
                                                POINTER(c_float * 9),
                                                POINTER(c_float * 3),
                                                POINTER(c_float * 9)]
    _TRICAL.TRICAL_estimate_get_ext.restype = None

    _TRICAL.TRICAL_measurement_calibrate.argtypes = [POINTER(_Instance),
                                                     POINTER(c_float * 3),
                                                     POINTER(c_float * 3)]
    _TRICAL.TRICAL_measurement_calibrate.restype = None


class Instance(object):
    def __init__(self, field_norm=1.0, measurement_noise=1e-6):
        """
        Create a new TRICAL instance with the supplied field norm (magnitude)
        and measurement noise.

        The field norm should be the expected magnitude of the field at the
        sensor location, and the measurement noise should be the standard
        deviation of the error in sensor readings (the error being presumed
        to be white Gaussian noise).
        """
        if not _TRICAL:
            _init()

        # Sanity-check the input parameters
        if field_norm <= 0.0:
            raise ValueError("Field norm must be > 0.0 (got %f)" % field_norm)

        if measurement_noise <= 0.0:
            raise ValueError("Measurement noise must be > 0.0 (got %f)" %
                             measurement_noise)

        # Initialize the internal (C) instance
        self._instance = _Instance()
        _TRICAL.TRICAL_init(self._instance)
        _TRICAL.TRICAL_norm_set(self._instance, field_norm)
        _TRICAL.TRICAL_noise_set(self._instance, measurement_noise)

        # Initialize the Python-accessible calibration estimate
        self.bias = (0.0, 0.0, 0.0)
        self.scale = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        self.measurement_count = 0

    def update(self, measurement):
        """
        Update the calibration estimate based on a new measurement.
        """
        if not measurement or len(measurement) != 3:
            raise ValueError("Measurement must be a sequence with 3 items")

        _TRICAL.TRICAL_estimate_update(self._instance,
                                       (c_float * 3)(*measurement))

        bias = (c_float * 3)()
        scale = (c_float * 9)()

        _TRICAL.TRICAL_estimate_get(self._instance, bias, scale)
        self.bias = tuple(bias[0:3])
        self.scale = tuple(scale[0:9])
        self.measurement_count = \
            _TRICAL.TRICAL_measurement_count_get(self._instance)

    def calibrate(self, measurement):
        """
        Given a measurement, return a calibrated measurement based on the
        current calibration estimate.
        """
        if not measurement or len(measurement) != 3:
            raise ValueError("Measurement must be a sequence with 3 items")

        calibrated_measurement = (c_float * 3)()
        _TRICAL.TRICAL_measurement_calibrate(self._instance,
                                             (c_float * 3)(*measurement),
                                             calibrated_measurement)

        return tuple(calibrated_measurement[0:3])


def _squared_norm(v):
    return v[0] * v[0] + v[1] * v[1] + v[2] * v[2]


def generate_html_viz(instance, samples):
    """
    Generate a WebGL visualisation of
    """
    import math
    import pkg_resources

    raw_data = samples
    iteratively_calibrated_data = []
    squared_magnitude = 0.0

    # Iterate over the samples and generate the iteratively calibrated data
    for measurement in samples:
        instance.update(measurement)
        calibrated_measurement = instance.calibrate(measurement)
        iteratively_calibrated_data.append(calibrated_measurement)

        # Update the maximum magnitude seen
        squared_magnitude = max(squared_magnitude, _squared_norm(measurement))
        squared_magnitude = max(squared_magnitude,
                                _squared_norm(calibrated_measurement))

    # And now insert them into the appropriate part of the HTML. Try the
    # pkg_resources way, but if the package hasn't actually been installed
    # then look for it in the same directory as this file
    try:
        html_path = pkg_resources.resource_filename("TRICAL", "viz.html")
        html = open(html_path, "rb").read()
    except IOError:
        html_path = os.path.join(os.path.dirname(__file__), "viz.html")
        html = open(html_path, "rb").read()

    # Insert the visualisation data -- output measurements as a single array
    # of points
    raw_points = ",\n".join(",".join("%.4f" % v for v in measurement)
                            for measurement in raw_data)
    calibrated_points = ",\n".join(",".join(
        "%.4f" % v for v in measurement) for measurement in calibrated_data)

    html = html.replace("{{raw}}", "[" + raw_points + "]")
    html = html.replace("{{calibrated}}", "[" + calibrated_points + "]")

    # TODO: maybe insert the calibration estimate as well?
    html = html.replace("{{magnitude}}",
                        "%.3f" % math.sqrt(squared_magnitude))
    html = html.replace("{{fieldNorm}}",
                        "%.3f" % instance._instance.field_norm)

    return html
