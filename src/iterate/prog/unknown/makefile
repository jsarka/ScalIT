#!/bin/bash

include ../../MakeMpi.inc

MOSB_LIBS  = ../../lib/libmosb.a  ../../lib/libosb.a

MMOD_OBJS = ../p_type/mosbtype.o  ../p_osb_pio/mosb.o ../p_osb_pio/mosbp.o 
PMOD_OBJS = ../p_type/mosbtype.o  ../p_osb_nopio/posb.o ../p_osb_nopio/posbp.o 

BINDIR=../../bin/osbexe

all:update mosb mosb1 posb

update:
	cp ../p_osb_pio/*.mod  ./
	cp ../p_osb_nopio/*.mod  ./

mosb: mosb.drv.o
	$(MPI_FC)  $(LDFLAGS) -o $(BINDIR)/m_iterate mosb.drv.o  $(MMOD_OBJS) $(MOSB_LIBS)
	rm -f mosb.drv.o
	echo 

mosb1: mosb1.drv.o
	$(MPI_FC)  $(LDFLAGS) -o $(BINDIR)/m1_iterate mosb1.drv.o  $(MMOD_OBJS) $(MOSB_LIBS)
	rm -f mosb.drv.o
	echo 

posb: posb.drv.o
	$(MPI_FC)  $(LDFLAGS) -o $(BINDIR)/p_iterate posb.drv.o $(PMOD_OBJS) $(MOSB_LIBS)
	rm -f posb.drv.o
	echo


clean:	
	rm -f *~ core 

cleanall:clean
	rm -f *.o *.mod *.MOD
	rm -f $(BINDIR)/m_iterate
	rm -f $(BINDIR)/m1_iterate
	rm -f $(BINDIR)/p_iterate
	rm -f $(BINDIR)/mosbEig 
	rm -f $(BINDIR)/mosbHX
	rm -f $(BINDIR)/mosbPX
	rm -f $(BINDIR)/mosbHij
	rm -f $(BINDIR)/mtest
	rm -f $(BINDIR)/mtest1
	rm -f $(BINDIR)/ptest

