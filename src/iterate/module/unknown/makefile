#!/bin/bash

include ../../Makefile.inc

OBJS = mosbtype.o ../../lib/libmosb.a  

all:mosbtype.o

mosbtype.o: mosbtype.f90      mosbtype.data.h  mosbtype.mosb.h \
            mosbtype.print.h  mosbtype.cal.h   

clean:	
	rm -f *.o *~ core *.mod

cleanall:clean

