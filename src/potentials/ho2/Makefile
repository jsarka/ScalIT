MOL := ho2
MOD := ja3

include $(PREFIX)/Makefile.setup

POT_OBJS := $(patsubst %.f90,%.o,$(wildcard *.f90)) $(patsubst %.f,%.o,$(wildcard *.f))

OSB_LIBS  := $(LIBPATH)/libham.$(LIBTYPE) $(LIBPATH)/libutil.a 
MOSB_LIBS := $(LIBPATH)/libham.$(LIBTYPE) $(LIBPATH)/libutil.a $(LIBPATH)/libutilmpi.a 

MOD_OBJS  := $(SRCPATH)/ham/$(MOD)/module/$(MOD).o  $(SRCPATH)/ham/$(MOD)/prog/$(MOD).drv.o
PMOD_OBJS := $(SRCPATH)/ham/$(MOD)/module/p$(MOD).o $(SRCPATH)/ham/$(MOD)/prog/p$(MOD).drv.o
MMOD_OBJS := $(SRCPATH)/ham/$(MOD)/module/m$(MOD).o $(SRCPATH)/ham/$(MOD)/prog/m$(MOD).drv.o

E_OBJS := $(POT_OBJS) $(SRCPATH)/ham/$(MOD)/index/index_AB2e.o  
O_OBJS := $(POT_OBJS) $(SRCPATH)/ham/$(MOD)/index/index_AB2o.o  
A_OBJS := $(POT_OBJS) $(SRCPATH)/ham/$(MOD)/index/index_ABC.o 

PSO_VLR := $(SRCPATH)/ham/$(MOD)/prog/psolr.drv.o
PSO_VBR := $(SRCPATH)/ham/$(MOD)/prog/psoBr.drv.o

.PHONY: all clean
all:                          \
$(BINPATH)/$(MOL)/$(MOL)vlr   \
$(BINPATH)/$(MOL)/$(MOL)vBR   \
$(BINPATH)/$(MOL)/$(MOL)      \
$(BINPATH)/$(MOL)/$(MOL)_e    \
$(BINPATH)/$(MOL)/$(MOL)_o    \
$(BINPATH)/$(MOL)/m$(MOL)     \
$(BINPATH)/$(MOL)/m$(MOL)_e   \
$(BINPATH)/$(MOL)/m$(MOL)_o   \
$(BINPATH)/$(MOL)/p$(MOL)     \
$(BINPATH)/$(MOL)/p$(MOL)_e   \
$(BINPATH)/$(MOL)/p$(MOL)_o   \
$(BINPATH)/$(MOL)/$(MOL)pot   \
$(BINPATH)/$(MOL)/$(MOL)fit

%.o: %.f90
	$(MPIFC) $(FFLAGS) -o $@ -c $^

%.o: %.f
	$(MPIFC) $(FFLAGS) -o $@ -c $^

clean:
	rm -f $(POT_OBJS)
	rm -f $(BINPATH)/$(MOL)/*

$(BINPATH)/$(MOL)/$(MOL)vlr: $(POT_OBJS)
	$(FC) -o $@ $(E_OBJS) $(PSO_VLR) $(OSB_LIBS) $(LDFLAGS)
	
$(BINPATH)/$(MOL)/$(MOL)vBR: $(POT_OBJS)
	$(FC)  -o $@ $(E_OBJS) $(PSO_VBR) $(OSB_LIBS) $(LDFLAGS)
	
$(BINPATH)/$(MOL)/$(MOL): $(POT_OBJS)
	$(FC)  -o $@ $(A_OBJS) $(MOD_OBJS) $(OSB_LIBS) $(LDFLAGS)
	
$(BINPATH)/$(MOL)/$(MOL)_e: $(POT_OBJS)
	$(FC)  -o $@ $(E_OBJS) $(MOD_OBJS) $(OSB_LIBS) $(LDFLAGS)
	
$(BINPATH)/$(MOL)/$(MOL)_o: $(POT_OBJS)
	$(FC)  -o $@ $(O_OBJS) $(MOD_OBJS) $(OSB_LIBS) $(LDFLAGS)
	
$(BINPATH)/$(MOL)/m$(MOL): $(POT_OBJS)
	$(MPIFC)  -o $@ $(A_OBJS) $(MMOD_OBJS) $(MOSB_LIBS) $(LDFLAGS)
	
$(BINPATH)/$(MOL)/m$(MOL)_e: $(POT_OBJS)
	$(MPIFC)  -o $@ $(E_OBJS) $(MMOD_OBJS) $(MOSB_LIBS) $(LDFLAGS)
	
$(BINPATH)/$(MOL)/m$(MOL)_o: $(POT_OBJS)
	$(MPIFC)  -o $@ $(O_OBJS) $(MMOD_OBJS) $(MOSB_LIBS) $(LDFLAGS)
	
$(BINPATH)/$(MOL)/p$(MOL): $(POT_OBJS)
	$(MPIFC)  -o $@ $(A_OBJS) $(PMOD_OBJS) $(MOSB_LIBS) $(LDFLAGS)
	
$(BINPATH)/$(MOL)/p$(MOL)_e: $(POT_OBJS)
	$(MPIFC) -o $@ $(E_OBJS) $(PMOD_OBJS) $(MOSB_LIBS) $(LDFLAGS)
	
$(BINPATH)/$(MOL)/p$(MOL)_o: $(POT_OBJS) 
	$(MPIFC) -o $@ $(O_OBJS) $(PMOD_OBJS) $(MOSB_LIBS) $(LDFLAGS)
	
$(BINPATH)/$(MOL)/$(MOL)pot: $(POT_OBJS)
	$(FC)  -o $@ $(E_OBJS) $(SRCPATH)/ham/$(MOD)/prog/pot.drv.o $(OSB_LIBS) $(LDFLAGS)
	
$(BINPATH)/$(MOL)/$(MOL)fit: $(POT_OBJS)
	$(FC)  -o $@ $(E_OBJS) $(SRCPATH)/ham/$(MOD)/prog/fitv.drv.o $(OSB_LIBS) $(LDFLAGS)
