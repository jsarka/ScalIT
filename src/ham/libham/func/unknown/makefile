#!/bin/bash
DOSB_LIB  = ../../lib/libdosb.a
FUNC_OBJS = getRoot.o  powell.o  linmin.o  simpleMin.o   \
       spline.o   splint.o  amoeba.o  simplex.o     

all:  $(FUNC_OBJS) 
	ar -rucs $(DOSB_LIB) $(FUNC_OBJS)
	rm -f $(FUNC_OBJS)	

test: powell simplex simplex1 

powell: powell.drv.o	
	$(FC) -o test_powell powell.drv.o $(DOSB_LIB)
	rm -f powell.drv.o

simple: simple.drv.o
	$(FC) -o test_simple  simple.drv.o $(DOSB_LIB)
	rm -f simple.drv.o

simplex: simplex.drv.o
	$(FC) -o test_simplex  simplex.drv.o $(DOSB_LIB)
	rm -f simplex.drv.o

simplex1: simplex1.drv.o
	$(FC) -o test_simplex1  simplex1.drv.o $(DOSB_LIB)
	rm -f simplex1.drv.o	

%.o:%.f90 
	$(FC) -c $<
#	$(FC) -c $*.f90 

clean:	
	rm -f *.o *~ core
	rm -f test_*

cleanall:clean
	rm -f *.mod *.MOD
