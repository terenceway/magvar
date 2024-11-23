#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "av_mag.h"

const char FILENAME[] = "WMM2020_TEST_VALUES.txt";

#define DEG2RAD (PI/180.0)
#define RAD2DEG (180.0/PI)
#define radians(x) ((x)*DEG2RAD)
#define degrees(x) ((x)*RAD2DEG)

int main() {
	FILE *fp;
	char buffer[1024];
	float max_error = 0.01;
	unsigned errors = 0, tests = 0;

	printf("TAP version 13\n");
	printf("1..1\n");

	fp = fopen(FILENAME, "r");
	if (fp == NULL) {
		perror(FILENAME);
		exit(EXIT_FAILURE);
	}

	while (fgets(buffer, sizeof(buffer), fp)) {
		float	date, height, latitude, longitude,
			declination, inclination,
			h, x, y, z, f, dD, dI, dH, dX, dY, dZ, dF;

		if (sscanf(buffer, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %g",
				&date, &height, &latitude, &longitude,
				&declination, &inclination,
				&h, &x, &y, &z, &f,
				&dD, &dI, &dH, &dX, &dY, &dZ, &dF) == 18) {

			float magvar = magnetic_variation(date,
					radians(latitude),
					radians(longitude),
					height * 1000.0);

			tests++;

			if (fabs(degrees(magvar) + declination) > 0.01) {
				max_error = fabs(degrees(magvar) + declination);
				errors ++;

				printf("not ok - expected %g, got %g\n",
				       declination, -degrees(magvar));
				break;				
			}
		}
	}

	fclose(fp);

	if (errors) {
		exit(EXIT_FAILURE);
		return 1;
	}

	printf("ok 1 - %d tests passed\n", tests);

	exit(EXIT_SUCCESS);
	return 0;
}
