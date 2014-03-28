#define _BSD_SOURCE

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <omp.h>
#include <mpi.h>

static void compute_threaded(double *K, const double *D, const double *I, int n, int m);

int rank;
int size;

int main(int argc, char **argv) {

	if (argc != 6) {
		fprintf(stderr, "Usage: parallel N M DTENSOR ITENSOR KTENSOR\n");
		return 1;
	}

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int n = atoi(argv[1]);
	int m = atoi(argv[2]);
	const char *dtensor_path = argv[3];
	const char *itensor_path = argv[4];
	const char *ktensor_path = argv[5];

	int m_local = (rank < size-1) ? m / size : m - (m / size) * rank;
	off_t ioffset = (off_t) (m / size) * rank * n * n;

	int dtensor_fd = open(dtensor_path, O_RDONLY);
	int itensor_fd = open(itensor_path, O_RDONLY);

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

	int ktensor_fd = -1;
	if (rank == 0) {
		ktensor_fd = open(ktensor_path, O_CREAT | O_RDWR, 0644);
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
	}

	double *dtensor = mmap(NULL, n * n * sizeof(double), 
	                       PROT_READ, MAP_PRIVATE, dtensor_fd, 0);
	double *itensor = mmap(NULL, m_local * n * n * sizeof(double),
	                       PROT_READ, MAP_PRIVATE, itensor_fd, ioffset * sizeof(double));

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

	double *ktensor = NULL;
	if (rank == 0) {
		ktensor = mmap(NULL, n * n * sizeof(double),
		               PROT_READ | PROT_WRITE, MAP_SHARED, ktensor_fd, 0);

		if (!ktensor) {
			fprintf(stderr, "error in mmap on ktensor: ");
			perror(NULL);
			fprintf(stderr, "\n");
			return 1;
		}
	}

	double *stensor = malloc(n * n * sizeof(double));
	compute_threaded(stensor, dtensor, itensor, n, m_local);

	munmap(dtensor, n * n * sizeof(double));
	munmap(itensor, m_local * n * n * sizeof(double));
	close(dtensor_fd);
	close(itensor_fd);

	// final summation
	MPI_Reduce(stensor, ktensor, n * n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	free(stensor);
	
	if (rank == 0) {
		munmap(ktensor, n * n * sizeof(double));
		close(ktensor_fd);
	}

	MPI_Finalize();
	return 0;
}

static void compute_threaded(double *K, const double *D, const double *I, int n, int m) {
	double *T = malloc(m * n * n * sizeof(double));
	double *S = NULL;
	int s;

	#pragma omp parallel shared(S, s)
	{
		int r = omp_get_thread_num();
		s = omp_get_num_threads();

		#pragma omp master
		{
			S = malloc(s * n * n * sizeof(double));
		}
		#pragma omp barrier

		// initialize S
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				S[r*n*n + i*n + j] = 0.0;
			}
		}

		#pragma omp for schedule(static)
		for (int A = 0; A < m; A++) {
	
			// T_ilA = D_kl I_ikA ; T[A][i][l] = D[l][k] * I[A][i][k]
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

			// S_ijr = T_ilA I_jlA ; K[i][j] = T[A][i][l] * I[A][j][l]
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					for (int l = 0; l < n; l++) {
						S[r*n*n + i*n + j] += T[A*n*n + i*n + l] * I[A*n*n + j*n + l];
					}
				}
			}
		}

		#pragma omp barrier
		#pragma omp master
		{
			free(T);
		}
	}
	
	// K_ij = S_ijr
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			double sum = 0.0;
			for (int r = 0; r < s; r++) {
				sum += S[r*n*n + i*n + j];
			}
			K[i*n + j] = sum;
		}
	}

	free(S);
}