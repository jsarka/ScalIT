
include $(PREFIX)/Makefile.setup

.PHONY: all clean

all: $(BINPATH)/iterate $(BINPATH)/m_iterate $(BINPATH)/m1_iterate $(BINPATH)/p_iterate $(BINPATH)/p2_iterate

clean:
	find . -name "*.o" -exec rm {} \;

%.o: %.f90
	$(MPIFC) $(FFLAGS) -o $@ -c $^

$(BINPATH)/iterate: $(SRCPATH)/iterate/module/a_osbtype.o $(SRCPATH)/iterate/module/osb.o $(SRCPATH)/iterate/module/osbw.o $(SRCPATH)/iterate/prog/osb.drv.o
	$(FC) -o $@ $^ $(LIBPATH)/libiter.a $(LIBPATH)/libutil.a $(LDFLAGS)

$(BINPATH)/m_iterate: $(SRCPATH)/iterate/module/a_mosbtype.o $(SRCPATH)/iterate/module/mosb.o $(SRCPATH)/iterate/module/mosbp.o $(SRCPATH)/iterate/prog/mosb.drv.o
	$(MPIFC) -o $@ $^ $(LIBPATH)/libitermpi.a $(LIBPATH)/libiter.a $(LIBPATH)/libutil.a $(LIBPATH)/libutilmpi.a $(LDFLAGS)

$(BINPATH)/m1_iterate: $(SRCPATH)/iterate/module/a_mosbtype.o $(SRCPATH)/iterate/module/mosb.o $(SRCPATH)/iterate/module/mosbp.o $(SRCPATH)/iterate/prog/mosb1.drv.o
	$(MPIFC) -o $@ $^ $(LIBPATH)/libitermpi.a $(LIBPATH)/libiter.a $(LIBPATH)/libutil.a $(LIBPATH)/libutilmpi.a $(LDFLAGS)

$(BINPATH)/p_iterate: $(SRCPATH)/iterate/module/a_mosbtype.o $(SRCPATH)/iterate/module/posb.o $(SRCPATH)/iterate/module/posbp.o $(SRCPATH)/iterate/prog/posb.drv.o
	$(MPIFC) -o $@ $^ $(LIBPATH)/libitermpi.a $(LIBPATH)/libiter.a $(LIBPATH)/libutil.a $(LIBPATH)/libutilmpi.a $(LDFLAGS)

$(BINPATH)/p2_iterate: $(SRCPATH)/iterate/module/a_mosbtype.o $(SRCPATH)/iterate/module/posb2.o $(SRCPATH)/iterate/module/posbp2.o $(SRCPATH)/iterate/prog/posb2.drv.o
	$(MPIFC) -o $@ $^ $(LIBPATH)/libitermpi.a $(LIBPATH)/libiter.a $(LIBPATH)/libutil.a $(LIBPATH)/libutilmpi.a $(LDFLAGS)

