#!/bin/bash

include ../../MakeMpi.inc

MOSB_LIB = ../../lib/libmosb.a

MORTH_OBJS = gsorth_mpi.o gsorth_cx_mpi.o gsorth_sx_mpi.o 

all: $(MORTH_OBJS)
	ar -rucs $(MOSB_LIB) $(MORTH_OBJS)
	rm -f $(MORTH_OBJS)

clean:	
	rm -f *~ core *.o

cleanall:clean

