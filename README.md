# The magvar module

This module provides the function

	`magnetic_variation(date, latitude, longitude, elevation)`

`date` must be a float - the Julian date.  Example: 2022.1 is in February.

`latitude` is in radians, positive for north, negative for south.

`longitude` is in radians, positive for east, negative for west.

`elevation` is in meters.

This returns what the geo-scientists would call 'declination' but what
aviators call 'magnetic variation'.

The result is in radians, and is negative for EAST "east is least" and
positive for WEST.

For example:

```
>>> from geomag import magnetic_variation
>>> from math import degrees, radians

>>> degrees(magnetic_variation(2022.05,
                               radians(25.7953611),
                               radians(-80.2901158),
                               3)
6.9474...
```

This shows the magnetic variation for the Miami Internation Airport (KMIA)
as 6.9 degrees WEST, which matches what the sectional chart for January 2022.



## Building from COF files

This distribution includes WMM2010.COF, WMM2015.COF, and WMM2020.COF.

These files were downloaded from the National Centers for Environmental
Information (NCEI), formerly known as the National Geophysical Data Center
(NGDC).

	[NCEI](https://www.ncei.noaa.gov/products/world-magnetic-model)

The file `wmm.c` contains static data compiled from these .COF files.

To generate the file from scratch, use the `compile-cof' script:

```
$ ./compile-cof WMM2010.COF WMM2015.COF WMM2020.COF >wmm.c
```

This incorporates all models, so you can get the magnetic variation
for any latitude and longitude from any date from January 1, 2010.

If you're tight on space, and only need the current magnetic model:

```
$ ./compile-cof WMM2020.COF >wmm.c
```


## Testing the Software

The `test` subdirectory has test values from the NCEI.

Test the Python:

```
$ cd test
$ python3 test_wmm2020_test_values.py
```

Test the C code:

```
$ cd test
$ make
$ ./test_geomag
$ make clean
```
