file:= dynamicOOP.x dynamic.x
CXXFLAGS:= -std=c++0x -pthread

all: ${file}

dynamicOOP.x: main.cxx Manager.o Node.o Node_master.o Node_slave.o
	mpicxx ${CXXFLAGS} $^ -o $@

%.o: %.cxx
	mpicxx ${CXXFLAGS} -c $< -o $@

%.x: %.cxx
	mpicxx ${CXXFLAGS} $< -o $@

.PHONY: test1 test2 clean

test1: ${file}
	./dynamicOOP.x 5 2

test2: ${file}
	./dynamic.x 3

clean:
	rm -rf *.x *.o set*
