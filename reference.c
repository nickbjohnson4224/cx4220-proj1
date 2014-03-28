#define _BSD_SOURCE

#include <stdlib.h>
#include <stdio.h>

#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

static void compute(double *K, const double *D, const double *I, int n);

int main(int argc, char **argv) {

	if (argc != 5) {
		fprintf(stderr, "Usage: reference N DTENSOR ITENSOR KTENSOR\n");
		return 1;
	}

	int n = atoi(argv[1]);
	const char *dtensor_path = argv[2];
	const char *itensor_path = argv[3];
	const char *ktensor_path = argv[4];

	int dtensor_fd = open(dtensor_path, O_RDONLY);
	int itensor_fd = open(itensor_path, O_RDONLY);
	int ktensor_fd = open(ktensor_path, O_CREAT | O_RDWR, 0644);

	if (dtensor_fd < 0) {
		fprintf(stderr, "error: %s: ", dtensor_path);
		perror(NULL);
		return 1;
	}

	if (itensor_fd < 0) {
		fprintf(stderr, "error: %s: ", itensor_path);
		perror(NULL);
		return 1;
	}

	if (ktensor_fd < 0) {
		fprintf(stderr, "error: %s: ", ktensor_path);
		perror(NULL);
		return 1;
	}

	if (ftruncate(ktensor_fd, n * n * n * 4 * sizeof(double)) < 0) {
		fprintf(stderr, "error: could not resize ktensor file: ");
		perror(NULL);
		return 1;
	}

	double *dtensor = mmap(NULL, n * n * sizeof(double), 
	                       PROT_READ, MAP_PRIVATE, dtensor_fd, 0);
	double *itensor = mmap(NULL, n * n * n * 4 * sizeof(double), 
	                       PROT_READ, MAP_PRIVATE, itensor_fd, 0);
	double *ktensor = mmap(NULL, n * n * n * 4 * sizeof(double),
	                       PROT_READ | PROT_WRITE, MAP_SHARED, ktensor_fd, 0);

	if (!dtensor) {
		fprintf(stderr, "error in mmap on dtensor: ");
		perror(NULL);
		fprintf(stderr, "\n");
		return 1;
	}

	if (!itensor) {
		fprintf(stderr, "error in mmap on itensor: ");
		perror(NULL);
		fprintf(stderr, "\n");
		return 1;
	}
	
	if (!ktensor) {
		fprintf(stderr, "error in mmap on ktensor: ");
		perror(NULL);
		fprintf(stderr, "\n");
		return 1;
	}

	compute(ktensor, dtensor, itensor, n);

	munmap(dtensor, n * n * sizeof(double));
	munmap(itensor, n * n * n * 4 * sizeof(double));
	munmap(ktensor, n * n * n * 4 * sizeof(double));

	close(dtensor_fd);
	close(itensor_fd);
	close(ktensor_fd);

	return 0;
}

static void compute(double *K, const double *D, const double *I, int n) {

	// T_ilA = D_kl I_ikA
	double *T = malloc(n * n * n * 4 * sizeof(double));
	for (int i = 0; i < n; i++) {
		for (int l = 0; l < n; l++) {
			for (int A = 0; A < 4*n; A++) {
				double sum = 0;
				for (int k = 0; k < n; k++) {
					// sum += D[k][l] * I[i][k][A]
					sum += D[k*n + l] * I[i*n*4*n + k*4*n + A];
				}
				// T[i][l][A] = sum
				T[i*n*4*n + l*4*n + A] = sum;
			}
		}
	}

	// K_ij = T_ilA I_jlA
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			double sum = 0;
			for (int l = 0; l < n; l++) {
				for (int A = 0; A < 4*n; A++) {
					// sum += T[i][l][A] * I[j][l][A]
					sum += T[i*n*4*n + l*4*n + A] * I[j*n*4*n + l*4*n + A];
				}
			}
			K[i*n + j] = sum;
		}
	}
}
