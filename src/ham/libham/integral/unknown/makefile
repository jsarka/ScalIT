#!/bin/bash

DOSB_LIB = ../../lib/libdosb.a

INT_OBJS = gauIChev.o  gauILege.o gauIHermit.o gauILagu.o gauIAssLagu.o \
           integralDVR.o integralFmk.o integralPjm.o

all:$(INT_OBJS)
	ar -rucs $(DOSB_LIB) $(INT_OBJS)
	rm -f $(INT_OBJS)

test:testGLege testGChev testGLagu testGALagu testGHerm

testGLege:gauLege.drv.o
	$(FC) -o test_GLege gauLege.drv.o $(DOSB_LIB) $(OSB_FLAGS)
	rm -f gauLege.drv.o
	echo

testGChev:gauChev.drv.o
	$(FC) -o test_GChev gauChev.drv.o $(DOSB_LIB) $(OSB_FLAGS)
	rm -f gauChev.drv.o
	echo 

testGLagu:gauLagu.drv.o 
	$(FC) -o test_GLagu gauLagu.drv.o $(DOSB_LIB) $(OSB_FLAGS)
	rm -f gauLagu.drv.o
	echo

testGALagu:gauAssLagu.drv.o 
	$(FC) -o test_GALagu gauAssLagu.drv.o $(DOSB_LIB) $(OSB_FLAGS)
	rm -f gauAssLagu.drv.o
	echo

testGHerm:gauHerm.drv.o 
	$(FC) -o test_GHerm gauHerm.drv.o $(DOSB_LIB) $(OSB_FLAGS)
	rm -f gauHerm.drv.o
	echo

%.o:%.f90
	$(FC) -c $<
#	$(FC) -c $*.f90 

clean:	
	rm -f *.o *~ core
	rm -f test_*

cleanall:clean
	rm -f *.mod *.MOD 

