#!/usr/bin/env python3

from math import degrees, radians

from magvar import magnetic_variation

files = 0

USAGE = "Usage: python3 test_values.py {TestValues text file}..."

def main(args):
    errors = 0

    if len(args) == 0:
        usage()

    print('TAP version 13')
    print('1..%d' % len(args))

    for fn in args:
        test_file(fn)

def test_file(fn):
    global files

    errors = count = 0

    files += 1

    with open(fn) as f:
        for line in f:
            a = line.split()

            if len(a) == 18:
                a = [float(x) for x in a]

                mv = magnetic_variation(a[0],          # Date
                                        radians(a[2]), # geodetic latitude
                                        radians(a[3]), # geodetic longitude
                                        a[1] * 1000)   # Height above ellipsoid

                error = abs(a[4] + degrees(mv))

                count += 1

                if error > 0.01:
                    print('not ok', files, '- expected', a[4], 'got', -degrees(mv))
                    errors += 1
                    break

    if not errors:
        print('ok 1 -', files, 'tests passed')

def usage():
    from sys import stderr, exit

    print(USAGE, file=stderr)
    exit(1)

if __name__ == '__main__':
    from sys import argv
    main(argv[1:])
