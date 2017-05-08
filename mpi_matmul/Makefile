CC:= mpicc
CXX:= mpicxx

FLAGS:= -fopenmp

all: matmul.x

matmul.x: driver.o matmul-mpi.o
	$(CXX) $(FLAGS) $^ -o $@

%.o: %.cc
	$(CC) $(FLAGS) -c $< -o $@

%.o: %.cxx
	$(CXX) $(FLAGS) -c $< -o $@

.PHONY: test clean

test: matmul.x prove.py
	mpirun -np 4 ./$< 16 8
	./$(word 2,$^) matrix_A.mat matrix_B.mat matrix_C.mat


clean:
	rm -rf *.o *.x *.mat
