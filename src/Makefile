CC = g++
FLAGS = -Wall -g -O2 -Wextra -fPIC

all: meshgen.exe Mesh blur

meshgen.exe:
	cd meshgen && make
	
FEM.o:
	$(CC) -c $(FLAGS) -o FEM.o ./Mesh/FEM.cpp
	
blur.o:
	$(CC) -c $(FLAGS) -o blur.o blur.cpp

blur: FEM.o blur.o
	$(CC) $(FLAGS) -o blur blur.o FEM.o


clean:
	cd meshgen && make clean && rm -rf mesh_output && cd .. && rm -rf blur *.o
