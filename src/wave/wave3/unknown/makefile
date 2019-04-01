#!/bin/bash

COBJS = ../../lib/libdosb.a  #dvr/getv.o ../index/index_pjm.o ../poly/yjm.o 

WAVE3_OBJS = wave3jb.o wave3hs.o hsjb.o

HJ_OBJS =  $(WAVE3_OBJS) $(COBJS)

BINDIR= ../../bin/wave3

all: $(WAVE3_OBJS)

test:all wave3

testall:test jbexch hsjb jbwf

wave3: wave3jb wave3jb_e wave3jb_o wave3hs  wave3hs_e  wave3hs_o   

WAVE3 = wave3.data wave3.init wave3.io

#--------------- Wave3JB Module --------------------------------
wave3jb.o:wave3jb.f90   wave3jb.init  wave3jb.io    wave3jb.wfp1  wave3jb.wfm1 \
          wave3jb.wfs1  wave3jb.wfp2  wave3jb.wfm2  wave3jb.wfs2  $(WAVE3)

wave3jb:$(WAVE3_OBJS) wave3jb.drv.o
	$(FC) -o $(BINDIR)/wave3jb  wave3jb.drv.o $(HJ_OBJS) ../index/index_ABC.o $(OSB_FLAGS)
	echo

wave3jb_e:$(WAVE3_OBJS) wave3jb.drv.o
	$(FC) -o $(BINDIR)/wave3jb_e wave3jb.drv.o $(HJ_OBJS) ../index/index_AB2e.o $(OSB_FLAGS)
	echo

wave3jb_o:$(WAVE3_OBJS) wave3jb.drv.o
	$(FC) -o $(BINDIR)/wave3jb_o wave3jb.drv.o $(HJ_OBJS) ../index/index_AB2o.o $(OSB_FLAGS)
	echo

#--------------------  Wave3HS Module  ---------------------------
wave3hs.o:wave3hs.f90   wave3hs.init  wave3hs.io    wave3hs.wfp1   wave3hs.wfm1 \
          wave3hs.wfs1  wave3hs.wfp2  wave3hs.wfm2  wave3hs.wfs2   $(WAVE3)

wave3hs:$(WAVE3_OBJS) wave3hs.drv.o
	$(FC) -o $(BINDIR)/wave3hs wave3hs.drv.o $(HJ_OBJS) ../index/index_ABC.o $(OSB_FLAGS)
	echo

wave3hs_e:$(WAVE3_OBJS) wave3hs.drv.o
	$(FC) -o $(BINDIR)/wave3hs_e wave3hs.drv.o $(HJ_OBJS) ../index/index_AB2e.o $(OSB_FLAGS)
	echo

wave3hs_o:$(WAVE3_OBJS) wave3hs.drv.o
	$(FC) -o $(BINDIR)/wave3hs_o wave3hs.drv.o $(HJ_OBJS) ../index/index_AB2o.o $(OSB_FLAGS)
	echo

#------------ Misc. Testings  --------------------------
jbexch:hsjb.o jbexch.drv.o
	$(FC) -o ./jbexch hsjb.o jbexch.drv.o  
	echo

hsjb:hsjb.o hsjb.drv.o
	$(FC) -o ./hsjb hsjb.o hsjb.drv.o
	echo

jbwf:jbwf.drv.o
	$(FC) -o ./jbwf jbwf.drv.o hsjb.o $(COBJS) ../index/index_AB2e.o $(OSB_FLAGS)
	echo

%.o:%.f90
	$(FC) -c $<
#       $(FC) -c $*.f90 

clean:
	rm -f *~ core
	rm -f ./jbexch
	rm -f ./hsjb
	rm -f ./jbwf

cleanall:clean
	rm -f *.o *.mod *.MOD 
	rm -f $(BINDIR)/wave3jb
	rm -f $(BINDIR)/wave3jb_e
	rm -f $(BINDIR)/wave3jb_o
	rm -f $(BINDIR)/wave3hs
	rm -f $(BINDIR)/wave3hs_e
	rm -f $(BINDIR)/wave3hs_o
