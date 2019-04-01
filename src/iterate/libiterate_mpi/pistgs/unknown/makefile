#!/bin/bash

include "../../MakeMpi.inc"

OSB_LIB = ../../lib/libosb.a
MOSB_LIB = ../../lib/libmosb.a

PISTAX = pist_mpi.o pistf_mpi.o pistHij_mpi.o pistEx_mpi.o    
PISTDX = pist_dx_mpi.o pistf_dx_mpi.o pistHij_dx_mpi.o pistEx_dx_mpi.o    
PISTCX = pist_cx_mpi.o pistf_cx_mpi.o pistHij_cx_mpi.o pistEx_cx_mpi.o     
PISTSX = pist_sx_mpi.o pistf_sx_mpi.o pistHij_sx_mpi.o pistEx_sx_mpi.o  

PIST = $(PISTAX) $(PISTCX) $(PISTDX) $(PISTSX)

all: $(PIST)
	ar -rucs $(MOSB_LIB) $(PIST)
	rm -f $(PIST)

test : test_pist test_pistdx test_pistcx test_pistsx   

test_pist:pist_mpi.drv.o
	$(MPI_FC) $(LDFLAGS)  -o test_pist pist_mpi.drv.o $(MOSB_LIB) $(OSB_LIB) $(OSB_FLAGS)
	rm -f pist_mpi.drv.o

test_pistcx:pist_cx_mpi.drv.o $(PISTCX) 
	$(MPI_FC)  $(LDFLAGS) -o test_pistcx pist_cx_mpi.drv.o $(MOSB_LIB) $(OSB_LIB) $(OSB_FLAGS)
	rm -f pist_cx_mpi.drv.o

test_pistsx:pist_sx_mpi.drv.o $(PISTSX) 
	$(MPI_FC)  $(LDFLAGS) -o test_pistsx pist_sx_mpi.drv.o $(MOSB_LIB) $(OSB_LIB) $(OSB_FLAGS)
	rm -f pist_sx_mpi.drv.o

test_pistdx:pist_dx_mpi.drv.o $(PISTDX)                      
	$(MPI_FC)  $(LDFLAGS)  -o test_pistdx pist_dx_mpi.drv.o $(MOSB_LIB) $(OSB_LIB) $(OSB_FLAGS)
	rm -f pist_dx_mpi.drv.o

%.o:%.f90
	$(FC) -c $<
#	$(FC) -c $*.f90 

clean:	
	rm -f *~ core *.o
	rm -f test_*

cleanall:clean
