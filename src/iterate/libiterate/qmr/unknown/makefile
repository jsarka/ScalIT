#!/bin/bash

include ../../Makefile.inc

OSB_LIB = ../../lib/libosb.a
QMR_OBJS = qmr.o qmrcx.o     

all: $(QMR_OBJS)
	ar -rcus $(OSB_LIB) $(QMR_OBJS)
	rm -f $(QMR_OBJS)

clean:	
	rm -f  *~ core *.o

cleanall:clean

