DATAGEN_OBJECTS := datagen.o
DATAGEN_SOURCES := datagen.c

REFERENCE_OBJECTS := reference.o
REFERENCE_SOURCES := reference.c

CFLAGS := -std=c99 -pedantic -Wall -Wextra -Werror
CFLAGS += -O3 -fomit-frame-pointer

TARGETS := datagen reference

all: $(TARGETS)

datagen: $(DATAGEN_OBJECTS)
	$(CC) $< -o $@

reference: $(REFERENCE_OBJECTS)
	$(CC) $< -o $@

clean:
	rm *.o $(TARGETS)
