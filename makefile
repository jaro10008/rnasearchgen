CC=g++
CFLAGS=-fopenmp -O3
LIBS=-lRNA -lgomp

.PHONY: all

all: RNASearchGen

RNASearchGen.o: RNASearchGen.cpp
	@$(CC) -c -o $@ $< $(CFLAGS)

RNASearchGen: RNASearchGen.o
	@$(CC) -o $@ $< $(LIBS)

.PHONY: clean

clean:
	rm -f RNASearchGen.o RNASearchGen
