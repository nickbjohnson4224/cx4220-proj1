#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>

int main(int argc, char **argv) {
	
	if (argc != 3) {
		fprintf(stderr, "Usage: ddiff TENSOR1 TENSOR2\n");
		return 1;
	}

	const char *tensor1_path = argv[1];
	const char *tensor2_path = argv[2];

	FILE *tensor1_file = fopen(tensor1_path, "rb");
	FILE *tensor2_file = fopen(tensor2_path, "rb");

	if (!tensor1_file) {
		fprintf(stderr, "error: %s: ", tensor1_path);
		perror(NULL);
		return 1;
	}

	if (!tensor2_file) {
		fprintf(stderr, "error: %s: ", tensor2_path);
		perror(NULL);
		return 1;
	}

	double epsilon = DBL_EPSILON;
	bool match = true;
	while (1) {
		double d1, d2;
		if (fread(&d1, sizeof(double), 1, tensor1_file) < 1) {
			if (fread(&d2, sizeof(double), 1, tensor2_file) < 1) {
				break;
			}
			match = false;
			break;
		}
		if (fread(&d2, sizeof(double), 1, tensor2_file) < 1) {
			match = false;
			break;
		}

		if (d1 == 0.0 && d2 == 0.0) {
			continue;
		}

		if (abs(d1 / d2 - 1) < epsilon) {
			continue;
		}

		printf("%f %f\n", d1, d2);
		match = false;
		break;
	}
	
	if (match) {
		printf("match\n");
	}
	else {
		printf("no match\n");
	}

	return 0;
}
