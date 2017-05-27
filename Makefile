CC= g++
CFLAGS = -std=c++11 -g -c -Wall -Wvla -Wextra -pthread -DNDEBUG
MFLAGS = -std=c++11 -Wall -Wvla -Wextra -O -pthread -DNDEBUG

#all targets
all: GenericMatrixDriver
	./GenericMatrixDriver

#the tar command
tar:
	tar cvf ex3.tar Matrix.hpp InitializationException.h \
			OperatorsException.h Makefile README

#executables

GenericMatrixDriver: GenericMatrixDriver.o Complex.o
	$(CC) $(MFLAGS) GenericMatrixDriver.o Complex.o -o GenericMatrixDriver

Matrix: Matrix.hpp InitializationException.h OperatorsException.h Complex.h
	$(CC) $(MFLAGS) Matrix.hpp -o Matrix.hpp.gch


#object files

GenericMatrixDriver.o: GenericMatrixDriver.cpp Complex.cpp
	$(CC) $(CFLAGS) GenericMatrixDriver.cpp -o GenericMatrixDriver.o

Complex.o: Complex.cpp
	$(CC) $(CFLAGS) Complex.cpp -o Complex.o

#other targets
clean:
	rm -f *.o ex3.tar Matrix.hpp.gch GenericMatrixDriver

	

