#!/bin/bash

OSB_LIB = ../lib/libosb.a
MOSB_LIB = ../lib/libmosb.a
MB_OBJS = index_node.o vector_mpi.o vector_cx_mpi.o   

all:$(MB_OBJS)
	ar -rucs $(MOSB_LIB) $(MB_OBJS)
	rm -f $(MB_OBJS)

test:test_vec test_veccx node nodelen

test_vec:vector_mpi.drv.o
	$(MPI_FC) -o test_vec vector_mpi.drv.o $(MOSB_LIB)
	rm -f vector_mpi.drv.o

test_veccx:vector_cx_mpi.drv.o
	$(MPI_FC) -o test_veccx vector_cx_mpi.drv.o $(OSB_LIB) $(MOSB_LIB)
	rm -f vector_cx_mpi.drv.o

node: node.drv.o
	$(MPI_FC) -o test_node node.drv.o $(MOSB_LIB)
	rm -f node.drv.o

nodelen:nodeLen.drv.o
	$(MPI_FC) -o test_len nodeLen.drv.o $(MOSB_LIB)
	rm -f nodeLen.drv.o

%.o:%.f90
	$(MPI_FC) -c $<
#	$(FC) -c $*.f90 

clean:	
	rm -f *~ core *.o
	rm -f test_*

cleanall:clean
	rm -f *.mod *.MOD

