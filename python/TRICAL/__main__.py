import os
import sys
import TRICAL

# Exit if incorrect arguments have been passed
if len(sys.argv) < 3:
    print "Usage: python -m TRICAL <field norm> <noise> [html]"
    sys.exit(1)

# Set up the instance
instance = TRICAL.Instance(field_norm=float(sys.argv[1]),
                           measurement_noise=float(sys.argv[2]))

if len(sys.argv) >= 4 and sys.argv[3] == 'html':
    # Run calibration and generate a WebGL visualisation of the raw and
    # calibrated measurements. Output the visualisation in HTML to stdout
    samples = []

    # Convert samples to floating-point values before passing them to
    # generate_html_viz
    for line in sys.stdin:
        measurement = map(float, line.strip("\n\r\t ").split(","))
        if len(measurement) != 3:
            continue
        samples.append(tuple(measurement))

    out = TRICAL.generate_html_viz(instance, samples)
    sys.stdout.write(out)
    sys.stdout.flush()
else:
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
