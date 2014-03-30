DATAGEN_OBJECTS := datagen.o
DATAGEN_SOURCES := datagen.c

DDIFF_OBJECTS := ddiff.o
DDIFF_SOURCES := ddiff.c

REFERENCE_OBJECTS := reference.o
REFERENCE_SOURCES := reference.c

PARALLEL_OBJECTS := parallel.o
PARALLEL_SOURCES := parallel.c

MPICC := mpicc

CFLAGS := -std=c99 -pedantic -Wall -Wextra
CFLAGS += -O3 -fomit-frame-pointer -march=native
CFLAGS += -fopenmp

TARGETS := datagen reference ddiff parallel-omp-mpi parallel-blas-mpi

all: $(TARGETS)

datagen: $(DATAGEN_OBJECTS)
	$(CC) $< -o $@

reference: $(REFERENCE_OBJECTS)
	$(CC) $< -o $@

threaded: $(THREADED_OBJECTS)
	$(CC) $< -o $@ -fopenmp

ddiff: $(DDIFF_OBJECTS)
	$(CC) $< -o $@

parallel-omp-mpi: $(PARALLEL_SOURCES)
	$(MPICC) $(CFLAGS) $< -o $@ -DUSE_OPENMP -DUSE_MPI

parallel-blas-mpi: $(PARALLEL_SOURCES)
	$(MPICC) $(CFLAGS) $< -o $@ -DUSE_BLAS -lopenblas -DUSE_MPI

#parallel-omp: $(PARALLEL_SOURCES)
#	$(CC) $(CFLAGS) $< -o $@ -DUSE_OPENMP

#parallel-blas: $(PARALLEL_SOURCES)
#	$(CC) $(CFLAGS) $< -o $@ -DUSE_BLAS -lopenblas

#%.o: %.c
#	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm *.o $(TARGETS)
