#!/bin/bash

include ../../Makefile.inc

OSB_LIB  = ../../lib/libosb.a
LAN_OBJS = lan.o lan_dx.o lanHij.o lanHij_dx.o lanEx.o lanEx_dx.o 

all:$(LAN_OBJS)
	ar -rucs $(OSB_LIB) $(LAN_OBJS) 
	rm -f $(LAN_OBJS)

clean:	
	rm -f *~ core *.o

cleanall:clean

