#!/bin/bash

include ../../MakeMpi.inc

MOSB_LIB = ../../lib/libmosb.a

PISTAX = pist_mpi.o pistf_mpi.o pistHij_mpi.o pistEx_mpi.o    
PISTDX = pist_dx_mpi.o pistf_dx_mpi.o pistHij_dx_mpi.o pistEx_dx_mpi.o    
PISTCX = pist_cx_mpi.o pistf_cx_mpi.o pistHij_cx_mpi.o pistEx_cx_mpi.o     
PISTSX = pist_sx_mpi.o pistf_sx_mpi.o pistHij_sx_mpi.o pistEx_sx_mpi.o  

PIST = $(PISTAX) $(PISTCX) $(PISTDX) $(PISTSX)

all: $(PIST)
	ar -rucs $(MOSB_LIB) $(PIST)
	rm -f $(PIST)

clean:	
	rm -f *~ core *.o

cleanall:clean

