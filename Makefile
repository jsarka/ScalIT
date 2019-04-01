# Get ScalIT root directory and check to see if any spaces are used
PREFIX := $(abspath $(dir $(word $(words $(MAKEFILE_LIST)),$(MAKEFILE_LIST))))
ifneq (1, $(words $(PREFIX)))
    $(error Fatal error: ScalIT cannot be built in a path that contains whitespace)
endif

# set output: these directories must exist already
INCPATH := $(PREFIX)/include
SRCPATH := $(PREFIX)/src
MODPATH := $(PREFIX)/module
LIBPATH := $(PREFIX)/lib
BINPATH := $(PREFIX)/bin

# get configuration information
include $(PREFIX)/Makefile.setup

# set compilation options
# eventually support fine-grain options, such as directory- or file-dependent compiler options
%.o: %.f90
	$(MPIFC) $(FFLAGS) -o $@ -c $^

#########################################	
########## overarching targets ##########

.PHONY: all clean pristine clean_lib_dir
all: binaries

pristine: clean pristine_binaries

clean: \
clean_lib_dir \
clean_util_lib \
clean_util_lib_mpi \
clean_presinc_module \
clean_ham_lib clean_ham_index clean_ham_drivers clean_ham_modules \
clean_iter_libs clean_iter_modules \
clean_binaries \
clean_wave3
	rm -f $(MODPATH)/*
	
clean_lib_dir:
	rm -f $(PREFIX)/lib/* 

#####################################	
########## utility library ##########

.PHONY: util_lib util_lib_mpi clean_util_lib clean_util_lib_mpi

util_lib: $(LIBPATH)/libutil.$(LIBTYPE)

LIBUTIL_OBJS := $(patsubst %.f90,%.o,$(shell find $(SRCPATH)/libutil/serial -name "*.f90"))
$(LIBPATH)/libutil.$(LIBTYPE): $(LIBUTIL_OBJS)
	if [ $(LIBTYPE) = a ]; then $(AR) -rcus $@ $^ ; else $(FC) -o $(LDFLAGS) $^ ; fi
	
util_lib_mpi: $(LIBPATH)/libutilmpi.$(LIBTYPE)
	
LIBUTIL_MPI_OBJS := $(patsubst %.f90,%.o,$(shell find $(SRCPATH)/libutil/parallel -name "*.f90"))
$(LIBPATH)/libutilmpi.$(LIBTYPE): $(LIBUTIL_MPI_OBJS)
	if [ $(LIBTYPE) = a ]; then $(AR) -rcus $@ $^ ; else $(FC) -o $(LDFLAGS) $^ ; fi
	
clean_util_lib:
	rm -f $(LIBUTIL_OBJS)
	
clean_util_lib_mpi:
	rm -f $(LIBUTIL_MPI_OBJS)

#####################################	
########## hamiltonians #############

.PHONY: ham_lib presinc_module ham_modules ham_drivers ham_index \
clean_ham_lib clean_presinc_module clean_ham_modules clean_ham_drivers clean_ham_index

# Hamiltonian library creation

ham_lib: presinc_module $(LIBPATH)/libham.$(LIBTYPE)

# Presinc module creation - Creates the PSOVBR data to be used in (JA3,JA4) modules

PRESINC_OBJS := $(patsubst %.f90,%.o,$(wildcard $(SRCPATH)/ham/presinc/*.f90))
presinc_module: $(PRESINC_OBJS)

# Specific files are added due to strange "find" methods.  TODO: fix this
LIBHAM_OBJS := $(SRCPATH)/ham/libham/coeff/gauntmod.o \
               $(SRCPATH)/ham/libham/coeff/threejmod.o \
               $(SRCPATH)/ham/libham/coeff/threej.o \
               $(SRCPATH)/ham/libham/ndpr/ndpr.o \
               $(patsubst %.f90,%.o,$(shell find $(SRCPATH)/ham/libham -name "*.f90"))
$(LIBPATH)/libham.$(LIBTYPE): $(PRESINC_OBJS) $(LIBHAM_OBJS)
	if [ $(LIBTYPE) = a ]; then $(AR) -rcus $@ $^ ; else $(FC) -o $(LDFLAGS) $^ ; fi
	
# (JA3,JA4) module creation - (3,4) atom systems described by Jacobian coords
	
JA3_HAM_MOD_OBJS := $(patsubst %.f90,%.o,$(wildcard $(SRCPATH)/ham/ja3/module/*.f90))
JA4_HAM_MOD_OBJS := $(patsubst %.f90,%.o,$(wildcard $(SRCPATH)/ham/ja4/module/*.f90))
ham_modules: $(JA3_HAM_MOD_OBJS) $(JA4_HAM_MOD_OBJS)

JA3_HAM_DRI_OBJS := $(patsubst %.f90,%.o,$(wildcard $(SRCPATH)/ham/ja3/prog/*.f90))
JA4_HAM_DRI_OBJS := $(patsubst %.f90,%.o,$(wildcard $(SRCPATH)/ham/ja4/prog/*.f90))
ham_drivers: $(JA3_HAM_DRI_OBJS) $(JA4_HAM_DRI_OBJS)

JA3_INDEX_OBJS := $(patsubst %.f90,%.o,$(wildcard $(SRCPATH)/ham/ja3/index/*.f90))
JA4_INDEX_OBJS := $(patsubst %.f90,%.o,$(wildcard $(SRCPATH)/ham/ja4/index/*.f90))
ham_index: $(JA3_INDEX_OBJS) $(JA4_INDEX_OBJS)
	
clean_ham_lib:
	rm -f $(LIBHAM_OBJS)
	
clean_presinc_module:
	rm -f $(PRESINC_OBJS)
	
clean_ham_modules:
	rm -f $(JA3_HAM_MOD_OBJS)
	rm -f $(JA4_HAM_MOD_OBJS)	
	
clean_ham_drivers:
	rm -f $(JA3_HAM_DRI_OBJS)
	rm -f $(JA4_HAM_DRI_OBJS)
	
clean_ham_index:
	rm -f $(JA3_INDEX_OBJS)
	rm -f $(JA4_INDEX_OBJS)
	
#####################################	
########## iterate ##################
.PHONY: iterate_lib iterate_lib_mpi iterate_module iterate_prog \
        clean_iter_libs clean_iter_modules
        
# library creation - serial and MPI
iterate_lib: $(LIBPATH)/libiter.$(LIBTYPE)
LIBITER_OBJS := $(patsubst %.f90,%.o,$(shell find $(SRCPATH)/iterate/libiterate -name "*.f90"))
$(LIBPATH)/libiter.$(LIBTYPE): $(LIBITER_OBJS)
	if [ $(LIBTYPE) = a ]; then $(AR) -rcus $@ $^ ; else $(FC) -o $(LDFLAGS) $^ ; fi

iterate_lib_mpi: $(LIBPATH)/libitermpi.$(LIBTYPE)
LIBITER_MPI_OBJS := $(patsubst %.f90,%.o,$(shell find $(SRCPATH)/iterate/libiterate_mpi -name "*.f90"))
$(LIBPATH)/libitermpi.$(LIBTYPE): $(LIBITER_MPI_OBJS)
	if [ $(LIBTYPE) = a ]; then $(AR) -rcus $@ $^ ; else $(FC) -o $(LDFLAGS) $^ ; fi
	
# module creation - serial, MPI with/without parallel IO (m/p)
ITER_MOD_OBJS := $(patsubst %.f90,%.o,$(wildcard $(SRCPATH)/iterate/module/*.f90))
iterate_module: $(ITER_MOD_OBJS)

ITER_DRI_OBJS := $(patsubst %.f90,%.o,$(wildcard $(SRCPATH)/iterate/prog/*.f90))
iterate_prog: $(ITER_DRI_OBJS)

clean_iter_libs:
	rm -f $(LIBITER_OBJS)
	rm -f $(LIBITER_MPI_OBJS)
	
clean_iter_modules:
	rm -f $(ITER_MOD_OBJS)
	rm -f $(ITER_DRI_OBJS) 
	
###########################################	
########## wavefunctions ##################

.PHONY: wave3 wave3_modules wave3_prog clean_wave3 create_wave_dirs

create_wave_dirs:
	@cd $(SRCPATH)/wave
	@for SYS in $(shell find $(SRCPATH)/wave/. -maxdepth 1 ! -path . -type d -printf '%f\n') ; do { \
	    if [ ! -d $(BINPATH)/$$SYS ]; then mkdir $(BINPATH)/$$SYS; fi; \
	}; done 

WAVE3_MOD_OBJS := $(patsubst %.f90,%.o,$(wildcard $(SRCPATH)/wave/wave3/module/hyperspherical/*.f90)) \
                  $(patsubst %.f90,%.o,$(wildcard $(SRCPATH)/wave/wave3/module/jacobi/*.f90))
wave3_modules: $(WAVE3_MOD_OBJS)

WAVE3_DRI_OBJS := $(patsubst %.f90,%.o,$(wildcard $(SRCPATH)/wave/wave3/prog/*.f90))
wave3_prog: $(WAVE3_DRI_OBJS)

wave3: create_wave_dirs wave3_modules wave3_prog
	$(MAKE) -C $(SRCPATH)/wave/wave3 \
	         PREFIX=$(PREFIX) SRCPATH=$(SRCPATH) BINPATH=$(BINPATH) \
            LIBPATH=$(LIBPATH) MODPATH=$(MODPATH) INCPATH=$(INCPATH) \
            all
            
clean_wave3:
	rm -f $(WAVE3_MOD_OBJS)
	rm -f $(WAVE3_DRI_OBJS)
	rm -f $(SRCPATH)/wave/wave3/*.o
	rm -rf $(BINPATH)/wave3

#############################################################	
########## binaries require everything to be built ##########
.PHONY: binaries ham iterate clean_binaries clean_iterate pristine_binaries

POTENTIALS := $(wildcard $(SRCPATH)/potentials/*)

binaries: create_dirs ham iterate wave3

create_dirs:
	@cd $(SRCPATH)/potentials
	@for SYS in $(shell find $(SRCPATH)/potentials/. -maxdepth 1 ! -path . -type d -printf '%f\n') ; do { \
	    if [ ! -d $(BINPATH)/$$SYS ]; then mkdir $(BINPATH)/$$SYS; fi; \
	}; done   

ham: util_lib util_lib_mpi ham_lib presinc_module ham_modules ham_drivers ham_index
	@for SYS in $(POTENTIALS); do { \
	    echo Building $$SYS system binaries ;\
	    $(MAKE) -C $$SYS -f Makefile \
                    PREFIX=$(PREFIX) SRCPATH=$(SRCPATH) BINPATH=$(BINPATH) \
                    LIBPATH=$(LIBPATH) MODPATH=$(MODPATH) INCPATH=$(INCPATH) \
                    all ;\
	    if [ $$? -ne 0 ]; then exit 1; fi;\
	}; done
	
iterate: util_lib util_lib_mpi iterate_lib iterate_lib_mpi iterate_module iterate_prog
	$(MAKE) -C $(SRCPATH)/iterate -f Makefile \
            PREFIX=$(PREFIX) SRCPATH=$(SRCPATH) BINPATH=$(BINPATH) \
            LIBPATH=$(LIBPATH) MODPATH=$(MODPATH) INCPATH=$(INCPATH) \
            all

clean_iterate:
	$(MAKE) -C $(SRCPATH)/iterate -f Makefile \
            PREFIX=$(PREFIX) SRCPATH=$(SRCPATH) BINPATH=$(BINPATH) \
            LIBPATH=$(LIBPATH) MODPATH=$(MODPATH) INCPATH=$(INCPATH) \
            clean
	
clean_binaries:
	@find $(SRCPATH)/potentials -name "*.o" -exec rm -f {} \;
	
pristine_binaries:
	@find $(BINPATH) -name .svn -prune -o -type f -exec rm -f {} \;
	rm -rf $(BINPATH)/*
