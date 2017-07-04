
file:= dynamic.x

all: ${file}

%.x: %.c
	mpicc $< -o $@

.PHONY: test

test: ${file}
	mpirun -np 1 $< 3
