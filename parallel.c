#define _BSD_SOURCE

#include <stdlib.h>
#include <malloc.h>
#include <stdint.h>
#include <stdio.h>

#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <fcntl.h>

#ifdef USE_OPENMP
#include <omp.h>
#endif

#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef USE_BLAS
#endif

#ifdef USE_BLAS
#define TRANSPOSE 't'
#define NO_TRANSPOSE 'n'
static void dgemm(double *C, double beta, const double *A, int An, int Am, char Atrans, 
                                          const double *B, int Bn, int Bm, char Btrans);
#endif
static double walltime(void);
static void compute_threaded(double *K, const double *D, const double *I, int n, int m);

int rank;
#ifdef USE_MPI
int size;
#endif

int main(int argc, char **argv) {

	if (argc != 6) {
		fprintf(stderr, "Usage: parallel N M DTENSOR ITENSOR KTENSOR\n");
		return 1;
	}

	#ifdef USE_MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	#else
	rank = 0;
	#endif

	int n = atoi(argv[1]);
	int m = atoi(argv[2]);
	const char *dtensor_path = argv[3];
	const char *itensor_path = argv[4];
	const char *ktensor_path = argv[5];

	#ifdef USE_MPI
	int m_local = (rank < size-1) ? m / size : m - (m / size) * rank;
	off_t ioffset = size * (m / size) * rank;
	#endif

	#ifdef USE_MPI
	printf("rank %d: m_local=%d ioffset=%d\n", rank, m_local, (int) ioffset);
	#endif

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

	double *dtensor = mmap(NULL, n * n * sizeof(double), 
	                       PROT_READ, MAP_PRIVATE | MAP_POPULATE, dtensor_fd, 0);
	
	#ifdef USE_MPI
	double *itensor = mmap(NULL, m_local * n * n * sizeof(double),
	                       PROT_READ, MAP_PRIVATE | MAP_POPULATE, itensor_fd, ioffset);
	#else
	double *itensor = mmap(NULL, m * n * n * sizeof(double),
                           PROT_READ, MAP_PRIVATE | MAP_POPULATE, itensor_fd, 0);
	#endif

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
	
	int ktensor_fd = -1;
	double *ktensor = NULL;
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
		
		ktensor = mmap(NULL, n * n * sizeof(double),
		               PROT_READ | PROT_WRITE, MAP_SHARED, ktensor_fd, 0);

		if (!ktensor) {
			fprintf(stderr, "error in mmap on ktensor: ");
			perror(NULL);
			fprintf(stderr, "\n");
			return 1;
		}
	}

	#ifdef USE_MPI
	double *stensor = malloc(n * n * sizeof(double));
	compute_threaded(stensor, dtensor, itensor, n, m_local);
	MPI_Reduce(stensor, ktensor, n * n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	#else
	compute_threaded(ktensor, dtensor, itensor, n, m);
	#endif

	munmap(dtensor, n * n * sizeof(double));
	munmap(itensor, m * n * n * sizeof(double));
	close(dtensor_fd);
	close(itensor_fd);

	if (rank == 0) {
		munmap(ktensor, n * n * sizeof(double));
		close(ktensor_fd);
	}

	#ifdef USE_MPI
	MPI_Finalize();
	#endif
	return 0;
}

#ifdef USE_BLAS
static void dgemm(double *C, double beta, 
                  const double *A, int An, int Am, char Atrans, 
                  const double *B, int Bn, int Bm, char Btrans) {
	int true_An = An;
	int true_Am = Am;
	int true_Bn = Bn;
	int true_Bm = Bm;

	if (Atrans == TRANSPOSE) {
		int t = Am;
		Am = An;
		An = t;
	}

	if (Btrans == TRANSPOSE) {
		int t = Bm;
		Bm = Bn;
		Bn = t;
	}

	if (true_Am != true_Bn) {
		fprintf(stderr, "error: attempt to multiply matrices %dx%d and %dx%d\n", 
			true_An, true_Am, true_Bn, true_Bm);
		exit(1);
	}

	double alpha = 1.0;
	dgemm_(&Atrans, &Btrans, &An, &Am, &Bm, &alpha, A, &true_An, B, &true_Bn, &beta, C, &An);
}
#endif

static double walltime(void) {
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return (double) (tp.tv_sec + tp.tv_usec * 1e-6);
}

static void compute_threaded(double *K, const double *D, const double *I, int n, int m) {
	#ifdef USE_MPI
	printf("rank %d: ", rank);
	#endif
	printf("compute_threaded %p %p %p %d %d\n", (void*) K, (void*) D, (void*) I, n, m);

	double tstart = walltime();

	#ifdef USE_OPENMP
	{
		// initialize K
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				K[i*n + j] = 0.0;
			}
		}

		#pragma omp parallel shared(K, D, I)
		{
			// initialize S matrix
			double S[n * n];
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					S[i*n + j] = 0.0;
				}
			}

			#pragma omp for schedule(static)
			for (int A = 0; A < m; A++) {
				double T[n];

				// T_ilA = D_kl I_ikA ; T[A][i][l] = D[l][k] * I[A][i][k]
				for (int i = 0; i < n; i++) {
					for (int l = 0; l < n; l++) {
						double sum = 0;
						for (int k = 0; k < n; k++) {
							// sum += D_kl * I_ikA ; sum += D[l][k] * I[A][i][k]
							sum += D[l*n + k] * I[A*n*n + i*n + k];
						}
						// T_ilA = sum ; T[A][i][l] = sum
						T[l] = sum;
					}

					// S_ijr = T_ilA I_jlA ; S_r[i][j] = T[A][i][l] * I[A][j][l]
					for (int j = 0; j < n; j++) {
						double sum = 0;
						for (int l = 0; l < n; l++) {
							sum += T[l] * I[A*n*n + j*n + l];
						}
						S[i*n + j] += sum;
					}
				}
			}

			// K_ij = S_ijr
			#pragma omp critical
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					K[i*n + j] += S[i*n + j];
				}
			}
		}
	}
	#endif//USE_OPENMP

	#ifdef USE_BLAS
	{
		double T[n * n];

		for (int A = 0; A < m; A++) {
			dgemm(T, 0.0, &I[A*n*n], n, n, NO_TRANSPOSE, D, n, n, TRANSPOSE);
			dgemm(K, (A==0) ? 0.0 : 1.0, T, n, n, NO_TRANSPOSE, &I[A*n*n], n, n, TRANSPOSE);
		}
	}
	#endif//USE_BLAS

	double tfinish = walltime();
	#ifdef USE_MPI
	printf("rank %d: ", rank);
	#endif
	printf("total compute time: %f seconds\n", tfinish - tstart);
}
