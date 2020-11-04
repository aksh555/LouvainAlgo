#!/bin/bash
XX= mpicxx
CC= g++

OPTFLAGS = -O3 -march=native -fopenmp -DPRINT_DIST_STATS -DDEBUG_PRINTF
SNTFLAGS = -std=c++11 -fopenmp -fsanitize=address -O1 -fno-omit-frame-pointer
CXXFLAGS = -std=c++11 -g $(OPTFLAGS)

ser= Serial
par1= OpenMP
par2= MPI
EXEC= community convert hierarchy
OBJ1= $(ser)/graph_binary.o $(ser)/community.o
OBJ2= $(ser)/graph.o
OBJS= $(par1)/utilities.o $(par1)/sorted-linked-list.o $(par1)/main.o $(par1)/input-handler.o $(par1)/execution-handler.o $(par1)/dynamic-weighted-graph.o $(par1)/community-exchange.o \
	  $(par1)/community-computation-weighted.o $(par1)/community-development.o
OBJ= $(par2)/main.o

all: serial openmp mpi
	@echo "Compile Success"

serial: community convert hierarchy
	@echo "Serial Success"

openmp: find-communities-OpenMP
	@echo "OpenMP Success"

mpi: detect
	@echo "MPI Success"

community : $(OBJ1) $(ser)/main_community.o
	$(CC) -o "community" $^

convert : $(OBJ2) $(ser)/main_convert.o
	$(CC) -o "convert" $^

hierarchy : $(ser)/main_hierarchy.o
	$(CC) -o "hierarchy" $^

%.o: community.cpp graph_binary.cpp graph.cpp main_community.cpp main_convert.cpp main_hierarchy.cpp
	$(CC) -o $@ -c $^

%.o: %.c
	gcc -fopenmp -O2 -w -c -o $@ $<


detect:
	$(XX) $(par2)/main.cpp $(CXXFLAGS) -c -o $(par2)/main.o
	$(XX) $(par2)/main.o $(OPTFLAGS) -o "detect"

find-communities-OpenMP: $(OBJS)
	gcc -fopenmp -w -o "find-communities-OpenMP" $^

clean:
	rm -f */*.o *.tree *level* */*.weights *.bin
	rm -f find-communities-OpenMP $(EXEC) detect
	@echo "Clean Success"
