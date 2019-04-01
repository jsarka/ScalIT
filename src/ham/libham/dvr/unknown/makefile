#!/bin/bash

DOSB_LIB = ../../lib/libdosb.a
DVR  = mtm.o mmt.o advr.o dvr.o dvr1.o dvr2.o dvrt.o dvrcx.o sdvr.o ddvr.o    #dvrCheck.o
SINC = sinc.o sinc_xh.o sinc1_dvr.o sinc2_dvr.o 

SINC_DVR = $(DVR) $(SINC)

OBJS = $(DVR) $(SINC) getv.o

all: $(OBJS) 
	ar -rucs $(DOSB_LIB) $(OBJS) 
	rm -f $(OBJS)

test: test_sinc1 test_sinc2 test_dvr1 test_dvr2

testall: test test_mm #test_gauss test_gausscx test_exp   

test_dvr1:dvr1.drv.o
	$(FC) -o test_dvr1 dvr1.drv.o $(DOSB_LIB) $(OSB_FLAGS)

test_dvr2:dvr2.drv.o
	$(FC) -o test_dvr2 dvr2.drv.o $(DOSB_LIB) $(OSB_FLAGS)

test_sinc1 : sinc1.drv.o
	$(FC) -o test_sinc1 sinc1.drv.o $(DOSB_LIB) $(OSB_FLAGS)

test_sinc2 : sinc2.drv.o
	$(FC) -o test_sinc2 sinc2.drv.o $(DOSB_LIB) $(OSB_FLAGS) 

test_gauss : gauss.drv.o 
	$(FC) -o test_gauss gauss.drv.o $(DOSB_LIB) $(OSB_FLAGS)

test_gausscx : gausscx.drv.o
	$(FC) -o test_gausscx gausscx.drv.o $(DOSB_LIB) $(OSB_FLAGS)

test_mm : mm.drv.o 
	$(FC) -o test_mm mm.drv.o $(DOSB_LIB) $(OSB_FLAGS)

test_exp: exp.drv.o 
	$(FC) -o test_exp  exp.drv.o  $(DOSB_LIB) $(OSB_FLAGS)

%.o:%.f90
	$(FC) -c $<
#	$(FC) -c $*.f90
 
clean:
	rm -f *.o *~ core
	rm -f test_*

cleanall:clean
          

