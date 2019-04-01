#!/bin/bash

include ../../Makefile.inc

OSB_LIB = ../../lib/libosb.a
AX_OBJS = pist.o    pistf.o    pistHij.o pistEx.o
DX_OBJS = pist_dx.o pistf_dx.o pistHij_dx.o pistEx_dx.o
CX_OBJS = pist_cx.o pistf_cx.o pistHij_cx.o pistEx_cx.o
SX_OBJS = pist_sx.o pistf_sx.o pistHij_sx.o pistEx_sx.o

OBJS = $(AX_OBJS) $(CX_OBJS) $(DX_OBJS) $(SX_OBJS)

all:$(OBJS)
	ar -rucs $(OSB_LIB) $(OBJS)
	rm -f $(OBJS)

clean:	
	rm -f *~ core *.o

cleanall:clean

