#include <stdlib.h>
#include <stdio.h>
#include <time.h>

static double random(void) {
	return (double) rand() / RAND_MAX;
}

int main(int argc, char **argv) {
	srand(time(NULL));

	if (argc != 4) {
		fprintf(stderr, "Usage: datagen N DTENSOR ITENSOR\n");
		return 1;
	}

	int n = atoi(argv[1]);
	const char *dtensor_path = argv[2];
	const char *itensor_path = argv[3];

	FILE *dtensor_file = fopen(dtensor_path, "wb");
	FILE *itensor_file = fopen(itensor_path, "wb");

	if (!dtensor_file) {
		fprintf(stderr, "error: %s:", dtensor_path);
		perror(NULL);
		return 1;
	}

	if  (!itensor_file) {
		fprintf(stderr, "error: %s:", itensor_path);
		perror(NULL);
		return 1;
	}

	for (int i = 0; i < n * n; i++) {
		double v = random();
		if (fwrite(&v, sizeof(double), 1, dtensor_file) < sizeof(double)) {
			fprintf(stderr, "error: writing to dtensor: ");
			perror(NULL);
			return 1;
		}
	}
	fclose(dtensor_file);

	for (int i = 0; i < n * n * n * 4; i++) {
		double v = random();
		if (fwrite(&v, sizeof(double), 1, itensor_file) < sizeof(double)) {
			fprintf(stderr, "error: writing to itensor: ");
			perror(NULL);
			return 1;
		}
	}
	fclose(itensor_file);

	return 0;
}
