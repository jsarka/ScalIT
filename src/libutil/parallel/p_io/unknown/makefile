#!/bin/bash

MOSB_LIB = ../lib/libmosb.a
MOBJS = mio.o mio_seq.o mio_diag.o mio_grid.o 
POBJS = pio_seq.o pio_diag.o pio_grid.o 
NOBJS = nio_seq.o nio_grid.o nio_diag.o

OBJS = $(MOBJS) $(POBJS) $(NOBJS)

all:$(OBJS)
	ar -rucs $(MOSB_LIB) $(OBJS)
	rm -f $(OBJS)
	echo

#test:test_mio test_nio test_pio
test:test_mdiag test_mgrid test_m0 test_nio

test_mio:test_mseq test_mdiag test_mgrid test_m0
test_nio:test_nseq test_ndiag test_ngrid
test_pio:test_pseq test_pdiag test_pgrid

#########  MPI parallel IO for very large data ###############
test_m0:mio.drv.o
	$(MPI_FC) -o test_mio mio.drv.o $(MOSB_LIB)
	rm -f mio.drv.o
	echo

test_mseq:mio_seq.drv.o 
	$(MPI_FC) -o test_mseq mio_seq.drv.o $(MOSB_LIB)
	rm -f mio_seq.drv.o
	echo

test_mdiag:mio_diag.drv.o 
	$(MPI_FC) -o test_mdiag mio_diag.drv.o $(MOSB_LIB)
	rm -f mio_diag.drv.o
	echo

test_mgrid:mio_grid.drv.o 
	$(MPI_FC) -o test_mgrid mio_grid.drv.o $(MOSB_LIB)
	rm -f mio_grid.drv.o
	echo

######## MPI parallel IO for small data size ################
test_nseq:nio_seq.drv.o
	$(MPI_FC) -o test_nseq nio_seq.drv.o $(MOSB_LIB)
	rm -f nio_seq.drv.o
	echo

test_ndiag:nio_diag.drv.o
	$(MPI_FC) -o test_ndiag nio_diag.drv.o $(MOSB_LIB)
	rm -f nio_diag.drv.o
	echo

test_ngrid:nio_grid.drv.o
	$(MPI_FC) -o test_ngrid nio_grid.drv.o $(MOSB_LIB)
	rm -f nio_grid.drv.o
	echo

######### MPI IO without parallel IO #########################
test_pseq:pio_seq.drv.o
	$(MPI_FC) -o test_pseq pio_seq.drv.o $(MOSB_LIB)
	rm -f pio_seq.drv.o
	echo

test_pdiag:pio_diag.drv.o
	$(MPI_FC) -o test_pdiag pio_diag.drv.o $(MOSB_LIB)
	rm -f nio_diag.drv.o
	echo

test_pgrid:pio_grid.drv.o
	$(MPI_FC) -o test_pgrid pio_grid.drv.o $(MOSB_LIB)
	rm -f pio_grid.drv.o
	echo

############################################################
test_dist:mio_dist.drv.o
	$(MPI_FC) -o test_dist  mio_dist.drv.o $(MOSB_LIB)
	rm -f mio_dist.drv.o
	echo

############################################################
%.o:%.f90
	$(FC) -c $<
#	$(FC) -c $*.f90 

clean:	
	rm -f *~ core *.o
	rm -f test_*

cleanall:clean
	rm -f *.mod *.MOD

