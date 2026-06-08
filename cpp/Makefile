CC=g++

all: DEP-GAP

# link
DEP-GAP: DEP-GAP.o rngstream.o
	$(CC) -fopenmp DEP-GAP.o rngstream.o -o DEP-GAP

# compile
DEP-GAP.o: DEP-GAP.cpp
	$(CC) -c -fopenmp DEP-GAP.cpp -w -std=c++11

rngstream.o: rngstream.cpp
	$(CC) -c -fopenmp -lgomp rngstream.cpp

# clean all files
clean:
	rm -f *.o;