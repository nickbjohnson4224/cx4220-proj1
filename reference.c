#define _BSD_SOURCE

#include <stdlib.h>
#include <stdio.h>

#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

static void compute(double *K, const double *D, const double *I, int n, int m);

int main(int argc, char **argv) {

	if (argc != 6) {
		fprintf(stderr, "Usage: reference N M DTENSOR ITENSOR KTENSOR\n");
		return 1;
	}

	int n = atoi(argv[1]);
	int m = atoi(argv[2]);
	const char *dtensor_path = argv[3];
	const char *itensor_path = argv[4];
	const char *ktensor_path = argv[5];

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

	if (ftruncate(ktensor_fd, n * n * sizeof(double)) < 0) {
		fprintf(stderr, "error: could not resize ktensor file: ");
		perror(NULL);
		return 1;
	}

	double *dtensor = mmap(NULL, n * n * sizeof(double), 
	                       PROT_READ, MAP_PRIVATE, dtensor_fd, 0);
	double *itensor = mmap(NULL, m * n * n * sizeof(double), 
	                       PROT_READ, MAP_PRIVATE, itensor_fd, 0);
	double *ktensor = mmap(NULL, n * n * sizeof(double),
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

	compute(ktensor, dtensor, itensor, n, m);

	munmap(dtensor, n * n * sizeof(double));
	munmap(itensor, m * n * n * sizeof(double));
	munmap(ktensor, n * n * sizeof(double));

	close(dtensor_fd);
	close(itensor_fd);
	close(ktensor_fd);

	return 0;
}

static void compute(double *K, const double *D, const double *I, int n, int m) {

	// T_ilA = D_kl I_ikA ; T[A][i][l] = D[l][k] * I[A][i][k]
	double *T = malloc(m * n * n * sizeof(double));
	for (int A = 0; A < m; A++) {
		for (int i = 0; i < n; i++) {
			for (int l = 0; l < n; l++) {
				double sum = 0;
				for (int k = 0; k < n; k++) {
					// sum += D_kl * I_ikA ; sum += D[l][k] * I[A][i][k]
					sum += D[l*n + k] * I[A*n*n + i*n + k];
				}
				// T_ilA = sum ; T[A][i][l] = sum
				T[A*n*n + i*n + l] = sum;
			}
		}
	}

	// K_ij = T_ilA I_jlA ; K[i][j] = T[A][i][l] * I[A][j][l]
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			double sum = 0.0;
			for (int A = 0; A < m; A++) {
				for (int l = 0; l < n; l++) {
					// sum += T_ilA * T_jlA ; sum += T[A][i][l] * I[A][j][l]
					sum += T[A*n*n + i*n + l] * I[A*n*n + j*n + l];
				}
			}
			// K_ij = sum ; K[i][j] = sum
			K[i*n + j] = sum;
		}
	}

	printf("K[0] = %f\n", K[0]);

	free(T);
}
