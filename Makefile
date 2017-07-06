
file:= dynamic.x
CXXFLAGS:= -std=c++0x -pthread

all: ${file}

%.x: %.c
	mpicc $< -o $@

%.x: %.cxx
	mpicxx ${CXXFLAGS} $< -o $@

.PHONY: test clean

test: ${file}
	mpirun -np 1 $< 3

clean:
	rm -rf *.x
