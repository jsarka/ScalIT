#!/bin/bash

include ../Makefile.inc

OSB_LIB   = ../lib/libosb.a
BASE_OBJS = osbtype.o mm.o diag.o vector.o vector_cx.o \
            matrix.o matrix_dx.o matrix_cx.o \
            matrix_sx.o util.o util_cx.o

all:$(BASE_OBJS)
	ar -rucs $(OSB_LIB) $(BASE_OBJS)
#	rm -f $(BASE_OBJS)

clean:	
	rm -f *~ core *.o
	rm -f test_mm test_util

cleanall:clean
	rm -f *.mod *.MOD


