#!/bin/bash

include ../../../Makefile.inc
      
OBJS = bj_diag.o bj_jacobi.o bj_vupdate.o ndprconv.o ndpr.o 

all: $(OBJS) 

ndpr.o:ndpr.f90 ndpr.data ndpr.io  ndpr.dpr ndpr.xyz ndpr.rtp ndpr.size 

.SUFFIXES: .f90 .o


clean:	
	rm -f test*  r3h0 res3 ndprconv *~ core

cleanall:clean
	rm -f *.o *.mod *.MOD 

