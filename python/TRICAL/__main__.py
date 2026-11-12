from __future__ import print_function

import sys

import TRICAL


def _samples(lines):
    """
    Parses comma-separated readings: 3 values per line (a measurement), or 6
    (a measurement followed by the expected field). Other lines are skipped.
    """
    for line in lines:
        try:
            values = [float(v) for v in line.strip("\n\r\t ").split(",")]
        except ValueError:
            continue

        if len(values) == 3:
            yield tuple(values), None
        elif len(values) == 6:
            yield tuple(values[0:3]), tuple(values[3:6])


def main(argv):
    # Exit if incorrect arguments have been passed
    if len(argv) < 3:
        print("Usage: python -m TRICAL <field norm> <noise> [html]",
              file=sys.stderr)
        return 1

    # Set up the instance
    instance = TRICAL.Instance(field_norm=float(argv[1]),
                               measurement_noise=float(argv[2]))

    if len(argv) >= 4 and argv[3] == "html":
        # Run calibration and generate a WebGL visualisation of the raw and
        # calibrated measurements. Output the visualisation in HTML to stdout
        samples = [(m, r) if r is not None else m
                   for m, r in _samples(sys.stdin)]

        sys.stdout.write(TRICAL.generate_html_viz(instance, samples))
        sys.stdout.flush()
    else:
        # Run calibration for each line in stdin, and print the calibrated
        # measurement to stdout
        for measurement, reference_field in _samples(sys.stdin):
            instance.update(measurement, reference_field)
            calibrated_measurement = instance.calibrate(measurement)

            out = ",".join(("%.7f" % m) for m in calibrated_measurement)
            sys.stdout.write(out + "\n")
            sys.stdout.flush()

        # Display a final calibration summary once done
        out = "################# CALIBRATION #################\n"
        out += " b = [%10.7f, %10.7f, %10.7f]\n" % instance.bias
        out += " D = [ [ %10.7f, %10.7f, %10.7f ]\n" % instance.scale[0:3]
        out += "       [ %10.7f, %10.7f, %10.7f ]\n" % instance.scale[3:6]
        out += "       [ %10.7f, %10.7f, %10.7f ] ]\n" % instance.scale[6:9]
        out += " measurements: %d, estimated noise: %.7f\n" % (
            instance.measurement_count, instance.noise_estimate)
        sys.stderr.write(out)
        sys.stderr.flush()

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
