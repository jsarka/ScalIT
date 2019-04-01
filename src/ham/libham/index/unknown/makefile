#!/bin/bash

DOSB_LIB = ../../lib/libdosb.a
IND_OBJS = index_pjm.o  index_CF.o   index_CFj.o   \
       index_CGmk.o index_CGj.o  index_CG.o    \
       index_3j0.o  index_gaunt.o index_node.o \
       index_nlm.o  index_jm.o index_node.o index_M1D.o \
       util.o       getm3.o    #index_cg0.o

A0_OBJS = index_ABC.o  index_AB2e.o index_AB2o.o  \
       index_ABCD.o index_A2B2_A1.o  index_A2B2_A1k.o  \

all: $(IND_OBJS) $(A0_OBJS) 
	ar -rucs $(DOSB_LIB) $(IND_OBJS)
	rm -f $(IND_OBJS)

test:pjm nlm jm1 jm2 cg cf gaunt threej0 node nodelen  #cg0

pjm:pjm.drv.o
	$(FC) -o testPjm pjm.drv.o $(DOSB_LIB)
	rm -f pjm.drv.o
	echo

nlm:nlm.drv.o
	$(FC) -o testNLM  nlm.drv.o $(DOSB_LIB)
	rm -f nlm.drv.o
	echo

jm1:jm1.drv.o
	$(FC) -o testJm1 jm1.drv.o $(DOSB_LIB)
	rm -f jm1.drv.o
	echo

jm2:jm2.drv.o
	$(FC) -o testJm2 jm2.drv.o $(DOSB_LIB)
	rm -f jm2.drv.o
	echo

cg:cg.drv.o
	$(FC) -o testCG cg.drv.o $(DOSB_LIB)
	rm -f cg.drv.o
	echo

cf:cf.drv.o
	$(FC) -o testCF cf.drv.o $(DOSB_LIB)
	rm -f cf.drv.o
	echo

gaunt:gaunt.drv.o
	$(FC) -o testga  gaunt.drv.o $(DOSB_LIB)
	rm -f gaunt.drv.o
	echo

threej0:threej0.drv.o
	$(FC) -o test3j0 threej0.drv.o  $(DOSB_LIB)
	rm -f threej0.drv.o
	echo

node:node.drv.o
	$(MPI_FC) -o testNode node.drv.o $(DOSB_LIB)
	rm -f node.drv.o
	echo

nodelen:nodeLen.drv.o
	$(MPI_FC) -o testlen nodeLen.drv.o $(DOSB_LIB)
	rm -f nodeLen.drv.o
	echo

m1d:m1d.drv.o
	$(MPI_FC) -o testm1d m1d.drv.o $(DOSB_LIB)
	rm -f m1d.drv.o
	echo

cg0:cg0.drv.o
	$(FC) -o testcg0 cg0.drv.o $(DOSB_LIB)
	rm -f cg0.drv.o
	echo


%.o:%.f90
	$(FC) -c $<
#	$(FC) -c $*.f90 

#%.o:%.f 
#	$(FC) -c $< 

clean:	
	rm -f test* *~ core

cleanall:clean
	rm -f *.o 

