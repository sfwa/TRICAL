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
readings on stdin (3 values per line, ending with \\n), and write iteratively
calibrated values on stdout in the same format.
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
    Loads the TRICAL library and sets up the cyptes interface.

    Called automatically the first time a Python instance is created.
    """

    global _TRICAL

    # TODO: search properly
    lib = os.path.join(os.path.dirname(__file__), "c", "libTRICAL.dylib")
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


if __name__ == "__main__":
    # Exit if incorrect arguments have been passed
    if len(sys.argv) != 3:
        print "Usage: python -m TRICAL.__init__ <field norm> <noise>"
        sys.exit(1)

    # Set up the instance
    instance = Instance(field_norm=float(sys.argv[1]),
                        measurement_noise=float(sys.argv[2]))

    # Run calibration for each line in stdin, and print the calibrated
    # measurement to stdout
    for line in sys.stdin:
        measurement = map(float, line.strip("\n\r\t ").split(","))
        if len(measurement) != 3:
            continue

        instance.update(measurement)
        calibrated_measurement = instance.calibrate(measurement)

        out = ",".join(("%.7f" % m) for m in calibrated_measurement)
        sys.stdout.write(out + "\n")
        sys.stdout.flush()

    # Display a final calibration summary once done
    out =  "################# CALIBRATION #################\n"
    out += " b = [%10.7f, %10.7f, %10.7f]\n" % instance.bias
    out += " D = [ [ %10.7f, %10.7f, %10.7f ]\n" % instance.scale[0:3]
    out += "       [ %10.7f, %10.7f, %10.7f ]\n" % instance.scale[3:6]
    out += "       [ %10.7f, %10.7f, %10.7f ] ]\n" % instance.scale[6:9]
    sys.stderr.write(out)
    sys.stderr.flush()
