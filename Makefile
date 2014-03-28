DATAGEN_OBJECTS := datagen.o
DATAGEN_SOURCES := datagen.c

DDIFF_OBJECTS := ddiff.o
DDIFF_SOURCES := ddiff.c

REFERENCE_OBJECTS := reference.o
REFERENCE_SOURCES := reference.c

THREADED_OBJECTS := threaded.o
THREADED_SOURCES := threaded.c

CFLAGS := -std=c99 -pedantic -Wall -Wextra -Werror
CFLAGS += -O3 -fomit-frame-pointer
CFLAGS += -fopenmp

TARGETS := datagen reference threaded ddiff

all: $(TARGETS)

datagen: $(DATAGEN_OBJECTS)
	$(CC) $< -o $@

reference: $(REFERENCE_OBJECTS)
	$(CC) $< -o $@

threaded: $(THREADED_OBJECTS)
	$(CC) $< -o $@ -fopenmp

ddiff: $(DDIFF_OBJECTS)
	$(CC) $< -o $@

clean:
	rm *.o $(TARGETS)
