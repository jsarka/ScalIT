#!/bin/bash

include ../../Makefile.inc

MOSB_LIB = ../../lib/libmosb.a
MUTIL_OBJS = seqdata.o    getsimpos.o   getpos.o     \
       sendnodes.o  recvnodes.o   sendindex.o    recvindex.o  \
       vindex.o     mvindex.o     mvxindex.o getga.o

all:update $(MUTIL_OBJS)
	ar -rucs $(MOSB_LIB) $(MUTIL_OBJS)
	rm -f $(MUTIL_OBJS)

update:
	cp ../p_type/mosbtype.mod ./

clean:	
	rm -f *.o *~ core 

cleanall:clean
	rm -f *.mod *.MOD 

