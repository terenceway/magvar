#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "av_mag.h"

const char FILENAME[] = "WMM2020_TEST_VALUES.txt";

#define DEG2RAD (PI/180.0)
#define RAD2DEG (180.0/PI)
#define radians(x) ((x)*DEG2RAD)
#define degrees(x) ((x)*RAD2DEG)

#define MAX_ERROR (0.01)

static unsigned errors = 0, files = 0;
static float max_error = 0.00;

static void test_file(const char *filename);

int main(int argc, char **argv) {
	int i;

	if (argc < 2) {
		fprintf(stderr, "Usage: test_magvar {TestValues text file}...\n");
		exit(EXIT_SUCCESS);
		return 0;
	}

	printf("TAP version 13\n");
	printf("1..%d\n", argc - 1);

	for (i = 1; i < argc; i++)
		test_file(argv[i]);

	exit(errors > 0 ? EXIT_FAILURE : EXIT_SUCCESS);
}

static void test_file(const char *filename) {
	FILE *fp = fopen(filename, "r");
	char buffer[1024];
	int result = 0;
	unsigned tests = 0;

	files ++;

	if (fp == NULL) {
		printf("not ok %u - ", files);
		perror(filename);
		return;
	}

	while (fgets(buffer, sizeof(buffer), fp)) {
		float	date, height, latitude, longitude,
			declination, inclination,
			h, x, y, z, f, dD, dI, dH, dX, dY, dZ, dF;

		float	error;

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

			error = fabs(degrees(magvar) + declination);

			if (error > max_error) max_error = error;

			if (error > MAX_ERROR) {
				errors ++;
				result ++;

				printf("not ok %u - expected %g, got %g\n",
				       files, declination, -degrees(magvar));
				break;				
			}
		}
	}

	if (result == 0)
		printf("ok %u - %u tests passed\n", files, tests);

	fclose(fp);
}
