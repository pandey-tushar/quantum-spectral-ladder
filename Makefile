CC ?= gcc
CFLAGS ?= -O3 -march=native -fopenmp
LIBS = -lflint -lmpfr -lgmp -lm

all: ladder2 ladder3

ladder2: solver/ladder2.c
	$(CC) $(CFLAGS) $< -o $@ $(LIBS)

ladder3: solver/ladder3.c
	$(CC) $(CFLAGS) $< -o $@ $(LIBS)

clean:
	rm -f ladder2 ladder3

.PHONY: all clean
