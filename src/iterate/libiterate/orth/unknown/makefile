#!/bin/bash

include ../../Makefile.inc
 
OSB_LIB   = ../../lib/libosb.a
ORTH_OBJS = gsorth.o gsorth_cx.o gsorth_sx.o     

all:$(ORTH_OBJS)
	ar -rucs $(OSB_LIB) $(ORTH_OBJS)
	rm -f $(ORTH_OBJS)

clean:	
	rm -f  *~ core *.o

cleanall:clean

