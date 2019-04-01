MOL := ne4
MOD := ja4

include $(PREFIX)/Makefile.setup

POT_OBJS := $(patsubst %.f90,%.o,$(wildcard *.f90)) $(patsubst %.f,%.o,$(wildcard *.f))

OSB_LIBS = $(LIBPATH)/libham.a $(LIBPATH)/libutil.a 
MOSB_LIBS= $(LIBPATH)/libham.a $(LIBPATH)/libutil.a $(LIBPATH)/libutilmpi.a 

MOD_OBJS  = $(SRCPATH)/ham/$(MOD)/module/$(MOD).o   $(SRCPATH)/ham/$(MOD)/prog/$(MOD).drv.o
MMOD_OBJS = $(SRCPATH)/ham/$(MOD)/module/m$(MOD).o  $(SRCPATH)/ham/$(MOD)/prog/m$(MOD).drv.o
PMOD_OBJS = $(SRCPATH)/ham/$(MOD)/module/p$(MOD).o  $(SRCPATH)/ham/$(MOD)/prog/p$(MOD).drv.o
 
AUU_OBJS = $(POT_OBJS) $(SRCPATH)/ham/$(MOD)/index/index_A2B2_Auu.o
AGG_OBJS = $(POT_OBJS) $(SRCPATH)/ham/$(MOD)/index/index_A2B2_Agg.o
AGU_OBJS = $(POT_OBJS) $(SRCPATH)/ham/$(MOD)/index/index_A2B2_Agu.o
AUG_OBJS = $(POT_OBJS) $(SRCPATH)/ham/$(MOD)/index/index_A2B2_Aug.o
A_OBJS   = $(POT_OBJS) $(SRCPATH)/ham/$(MOD)/index/index_ABCD.o 

PSO_VLR1 = $(SRCPATH)/ham/$(MOD)/prog/psolr1.drv.o
PSO_VLR2 = $(SRCPATH)/ham/$(MOD)/prog/psolr2.drv.o
PSO_VBR  = $(SRCPATH)/ham/$(MOD)/prog/psoBr.drv.o

.PHONY: all clean

all: $(POT_OBJS)                      \
     $(BINPATH)/$(MOL)/$(MOL)vlr1     \
     $(BINPATH)/$(MOL)/$(MOL)vlr2     \
     $(BINPATH)/$(MOL)/$(MOL)vBR      \
     $(BINPATH)/$(MOL)/$(MOL)         \
     $(BINPATH)/$(MOL)/$(MOL)_Agg     \
     $(BINPATH)/$(MOL)/$(MOL)_Auu     \
     $(BINPATH)/$(MOL)/$(MOL)_Agu     \
     $(BINPATH)/$(MOL)/$(MOL)_Aug     \
     $(BINPATH)/$(MOL)/p$(MOL)        \
     $(BINPATH)/$(MOL)/p$(MOL)_Agg    \
     $(BINPATH)/$(MOL)/p$(MOL)_Auu    \
     $(BINPATH)/$(MOL)/p$(MOL)_Agu    \
     $(BINPATH)/$(MOL)/p$(MOL)_Aug    \
     $(BINPATH)/$(MOL)/m$(MOL)        \
     $(BINPATH)/$(MOL)/m$(MOL)_Agg    \
     $(BINPATH)/$(MOL)/m$(MOL)_Auu    \
     $(BINPATH)/$(MOL)/m$(MOL)_Agu    \
     $(BINPATH)/$(MOL)/m$(MOL)_Aug    

                 
clean:
	rm -f $(POT_OBJS)
	rm -f $(BINPATH)/$(MOL)/*

%.o: %.f90
	$(MPIFC) -o $@ $(FFLAGS) -c $^

%.o: %.f
	$(MPIFC) -o $@ $(FFLAGS) -c $^

$(BINPATH)/$(MOL)/$(MOL)vlr1: $(POT_OBJS)
	$(FC) -o $@ $(A_OBJS) $(PSO_VLR1) $(OSB_LIBS) $(LDFLAGS)
	
$(BINPATH)/$(MOL)/$(MOL)vlr2: $(POT_OBJS)
	$(FC) -o $@ $(A_OBJS) $(PSO_VLR2) $(OSB_LIBS) $(LDFLAGS)

$(BINPATH)/$(MOL)/$(MOL)vBR: $(POT_OBJS)
	$(FC) -o $@ $(A_OBJS) $(PSO_VBR) $(OSB_LIBS) $(LDFLAGS)
	
$(BINPATH)/$(MOL)/$(MOL): $(POT_OBJS)
	$(FC) -o $@ $(A_OBJS) $(MOD_OBJS) $(OSB_LIBS) $(LDFLAGS)

$(BINPATH)/$(MOL)/$(MOL)_Agg: $(POT_OBJS)
	$(FC) -o $@ $(AGG_OBJS) $(MOD_OBJS) $(OSB_LIBS) $(LDFLAGS)

$(BINPATH)/$(MOL)/$(MOL)_Auu: $(POT_OBJS)
	$(FC) -o $@ $(AUU_OBJS) $(MOD_OBJS) $(OSB_LIBS) $(LDFLAGS)

$(BINPATH)/$(MOL)/$(MOL)_Agu: $(POT_OBJS)
	$(FC) -o $@ $(AGU_OBJS) $(MOD_OBJS) $(OSB_LIBS) $(LDFLAGS)

$(BINPATH)/$(MOL)/$(MOL)_Aug: $(POT_OBJS)
	$(FC) -o $@ $(AUG_OBJS) $(MOD_OBJS) $(OSB_LIBS) $(LDFLAGS)

$(BINPATH)/$(MOL)/p$(MOL): $(POT_OBJS)
	$(MPIFC) -o $@ $(A_OBJS) $(PMOD_OBJS) $(MOSB_LIBS) $(LDFLAGS)

$(BINPATH)/$(MOL)/p$(MOL)_Agg: $(POT_OBJS)
	$(MPIFC) -o $@ $(AGG_OBJS) $(PMOD_OBJS) $(MOSB_LIBS) $(LDFLAGS)

$(BINPATH)/$(MOL)/p$(MOL)_Auu: $(POT_OBJS)
	$(MPIFC) -o $@ $(AUU_OBJS) $(PMOD_OBJS) $(MOSB_LIBS)  $(LDFLAGS)

$(BINPATH)/$(MOL)/p$(MOL)_Agu: $(POT_OBJS)
	$(MPIFC) -o $@ $(AGU_OBJS) $(PMOD_OBJS) $(MOSB_LIBS)  $(LDFLAGS)

$(BINPATH)/$(MOL)/p$(MOL)_Aug: $(POT_OBJS)
	$(MPIFC) -o $@ $(AUG_OBJS) $(PMOD_OBJS) $(MOSB_LIBS)  $(LDFLAGS)

$(BINPATH)/$(MOL)/m$(MOL): $(POT_OBJS)
	$(MPIFC) -o $@ $(A_OBJS) $(MMOD_OBJS) $(MOSB_LIBS)  $(LDFLAGS)

$(BINPATH)/$(MOL)/m$(MOL)_Agg: $(POT_OBJS)
	$(MPIFC) -o $@ $(AGG_OBJS) $(MMOD_OBJS) $(MOSB_LIBS) $(LDFLAGS)

$(BINPATH)/$(MOL)/m$(MOL)_Auu: $(POT_OBJS)
	$(MPIFC) -o $@ $(AUU_OBJS) $(MMOD_OBJS) $(MOSB_LIBS) $(LDFLAGS)

$(BINPATH)/$(MOL)/m$(MOL)_Agu: $(POT_OBJS)
	$(MPIFC) -o $@ $(AGU_OBJS) $(MMOD_OBJS) $(MOSB_LIBS) $(LDFLAGS)

$(BINPATH)/$(MOL)/m$(MOL)_Aug: $(POT_OBJS)
	$(MPIFC) -o $@ $(AUG_OBJS) $(MMOD_OBJS) $(MOSB_LIBS) $(LDFLAGS)

