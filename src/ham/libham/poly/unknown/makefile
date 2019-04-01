#!/bin/bash

DOSB_LIB  = ../../lib/libdosb.a
POLY_OBJS = chev.o pj.o yj.o pjm.o yjm.o poly.o  lagu.o  \
     trigo.o trigonom.o gauChev.o  gauLege.o gauHermit.o gauLagu.o

all:$(POLY_OBJS)
	ar -rucs $(DOSB_LIB) $(POLY_OBJS)
	rm -f $(POLY_OBJS)

test:testch testlagu testpj testyj testpjm testyjm testpy \
     testtri   testGLege

testch:chev.drv.o
	$(FC) -o testch chev.drv.o $(DOSB_LIB)
	rm -f .drv.o
	echo

testtri:trigo.drv.o 
	$(FC) -o testtri trigo.drv.o $(DOSB_LIB)
	rm -f trigo.drv.o
	echo

testlagu:lagu.drv.o 
	$(FC) -o testla lagu.drv.o $(DOSB_LIB) $(OSB_FLAGS)
	rm -f lagu.drv.o
	echo

testpj:pj.drv.o
	$(FC) -o testpj pj.drv.o $(DOSB_LIB)
	rm -f pj.drv.o
	echo

testyj:yj.drv.o 
	$(FC) -o testyj yj.drv.o $(DOSB_LIB)
	rm -f yj.drv.o
	echo

testpjm:pjm.drv.o 
	$(FC) -o testpjm pjm.drv.o $(DOSB_LIB)
	rm -f pjm.drv.o
	echo

testyjm:yjm.drv.o 
	$(FC) -o testyjm yjm.drv.o $(DOSB_LIB)
	rm -f yjm.drv.o
	echo

testpy:pyjm.drv.o 
	$(FC) -o testpy pyjm.drv.o $(DOSB_LIB)
	rm -f pyjm.drv.o
	echo

testGLege:gauLege.drv.o 
	$(FC) -o testGLege gauLege.drv.o $(DOSB_LIB) $(OSB_FLAGS) 
	rm -f gauLege.drv.o
	echo

%.o:%.f90
	$(FC) -c $<
#	$(FC) -c $*.f90 

clean:	
	rm -f *.o *~ core
	rm -f test*

cleanall:clean
	rm -f *.mod *.MOD 

