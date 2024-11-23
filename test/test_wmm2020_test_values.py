#!/usr/bin/env python3

from math import degrees, radians

from magvar import magnetic_variation

def main():
    errors = 0
    count = 0

    print('TAP version 13')
    print('1..1')

    with open('WMM2020_TEST_VALUES.txt') as f:
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
                    print('not ok - expected', a[4], 'got', -degrees(mv))
                    errors += 1
                    break



    if errors:
        from sys import exit
        exit(1)

    print('ok 1 -', count, 'tests passed')

if __name__ == '__main__':
    main()
