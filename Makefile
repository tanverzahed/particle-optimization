
CC = icpc
MPCC = mpicxx
OPENMP = -qopenmp
CFLAGS = -O3
LIBS =
MPILIBS =

TARGETS = serial openmp mpi autograder

all:	$(TARGETS)

serial: serial.o common.o
	$(CC) -o $@ $(LIBS) serial.o common.o
autograder: autograder.o common.o
	$(CC) -o $@ $(LIBS) autograder.o common.o
openmp: openmp.o common.o
	$(CC) -o $@ $(LIBS) $(OPENMP) openmp.o common.o
mpi: mpi.o common_mpi.o
	$(MPCC) -o $@ $(LIBS) $(MPILIBS) mpi.o common.o

autograder.o: autograder.cpp common.h
	$(CC) -c $(CFLAGS) autograder.cpp
openmp.o: openmp.cpp common.h
	$(CC) -c $(OPENMP) $(CFLAGS) openmp.cpp
serial.o: serial.cpp common.h
	$(CC) -c $(CFLAGS) serial.cpp
mpi.o: mpi.cpp common.h
	$(MPCC) -c $(CFLAGS) mpi.cpp -lifcore
common.o: common.cpp common.h
	$(CC) -c $(CFLAGS) common.cpp
common_mpi.o: common.cpp common.h
	$(MPCC) -c $(CFLAGS) common.cpp
clean:
	rm -f *.o $(TARGETS) *.stdout *.txt
