#!/bin/bash
include ../../Makefile.inc

MOSB_LIB  = ../../lib/libmosb.a
MQMR_OBJS = qmr_mpi.o qmrcx_mpi.o     

all: $(MQMR_OBJS)
	ar -rucs $(MOSB_LIB) $(MQMR_OBJS)
	rm -f $(MQMR_OBJS)

clean:	
	rm -f *~ core *.o

cleanall:clean

