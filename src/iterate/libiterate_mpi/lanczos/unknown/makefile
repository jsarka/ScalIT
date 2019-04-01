#!/bin/bash

include ../../MakeMpi.inc

MOSB_LIB = ../../lib/libmosb.a
OTHERS=../../p_base/vector_mpi.o ../../obase/util.o ../orth/gsorth_mpi.o
CX_OTHERS=../../p_base/vector_cx_mpi.o ../../obase/util.o ../../obase/util_cx.o ../orth/gsorth_cx_mpi.o

MLAN_OBJS = lan_mpi.o lan_cx_mpi.o lan_dx_mpi.o \
       lanHij_mpi.o lanHij_dx_mpi.o lanEx_mpi.o lanEx_dx_mpi.o 

all:$(MLAN_OBJS)
	ar -rucs $(MOSB_LIB) $(MLAN_OBJS)
#	rm -f $(MLAN_OBJS)

test : all test_lan     test_landx   

test_lan:lan_mpi.o  lanHij_mpi.o lan_mpi.drv.o
	$(MPI_FC) $(FCFLAGS) $(OTHERS)  -o test_lan lan_mpi.o lanHij_mpi.o lanEx_mpi.o lan_mpi.drv.o 

test_landx:lan_dx_mpi.o lanHij_dx_mpi.o lan_dx_mpi.drv.o 
	$(MPI_FC) $(FCFLAGS) $(CX_OTHERS) -o test_landx lan_dx_mpi.o lanHij_dx_mpi.o lanEx_dx_mpi.o lan_dx_mpi.drv.o 

clean:	
	rm -f *~ core *.o

cleanall:clean
	rm -f *.o 

