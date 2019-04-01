#!/bin/bash

DOSB_LIB = ../../lib/libdosb.a
#MOD_DIR  = ../../mod

OBJS = presinc.o

all:update test

update: $(OBJS)  
	ar -rcus $(DOSB_LIB) $(OBJS) 
	rm -f $(OBJS)
	echo
# mv presinc.mod $(MOD_DIR)

test: pist show print diffDat  # test_psinc

psinc: presinc.drv.o
	$(FC) -o test_psinc  presinc.o presinc.drv.o $(DOSB_LIB) $(OSB_FLAGS)
	rm -f presinc.drv.o
	echo

pist: pist.drv.o 
	$(FC) -o ../../bin/pist  pist.drv.o $(OSB_FLAGS)
	rm -f pist.drv.o
	echo

show: showData.drv.o 
	$(FC) -o ../../bin/showData showData.drv.o	
	rm -f showData.drv.o
	echo

print: printData.drv.o 
	$(FC) -o ../../bin/printData printData.drv.o	
	rm -f printData.drv.o
	echo

diffDat: diffData.drv.o 
	$(FC) -o ../../bin/diffData diffData.drv.o	
	rm -f diffData.drv.o
	echo

%.o:%.f90
	$(FC) -c $<
#	$(FC) -c $*.f90 

clean:	
	rm -f test* *~ core
	rm -f ../../bin/pist ../../bin/showData 
	rm -f ../../bin/printData ../../bin/diffData

cleanall:clean
	rm -f *.o *.mod *.MOD 

