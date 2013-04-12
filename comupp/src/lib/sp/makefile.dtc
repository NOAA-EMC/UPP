#--------------------------------------------------------------------
# This makefile is for the sp library
#
# The makefile expect two variables set from the command line 
#--------------------------------------------------------------------

SHELL = /bin/sh
LIB   = libsp.a

#
# build configuration determined before compile
include ../../../configure.upp

#
# extra flags
EXTRA_FFLAGS  = -c $(PROMOTION) -I$(INCMOD)
EXTRA_CFLAGS  = -c
EXTRA_ARFLAGS =

#
# What files are we looking for - 
# Note the new source for files w/ .F requiring preprocessing
OBJS_F   = bll2ps.o spanaly.o spdz2uv.o speps.o spfft1.o       \
           spffte.o spfft.o spfftpt.o spgradq.o spgradx.o spgrady.o    \
           splaplac.o splegend.o sppad.o spsynth.o sptezd.o            \
           sptez.o sptezmd.o sptezm.o sptezmv.o sptezv.o sptgpmd.o     \
           sptgpm.o sptgpmv.o sptgpsd.o sptgps.o sptgpsv.o sptgptd.o   \
           sptgpt.o sptgptsd.o sptgptvd.o sptgptv.o sptrand.o sptran.o \
           sptranf0.o sptranf1.o sptranf.o sptranfv.o sptranv.o        \
           sptrund.o sptrun.o sptrung.o sptrungv.o sptrunl.o sptrunm.o \
           sptrunmv.o sptruns.o sptrunsv.o sptrunv.o spuv2dz.o         \
           spvar.o spwget.o
OBJS1_F  = splat.o ncpus.o lapack_gen.o fftpack.o
OBJS_F77 = 
OBJS_C   = 
OBJS     = $(OBJS_F) $(OBJS1_F) $(OBJS_F77) $(OBJS_C)

SRCS_F   = $(OBJS_F:.o=.f)
SRCS1_F  = $(OBJS1_F:.o=.F)
SRCS_F77 = $(OBJS_F77:.o=.f)
SRCS_C   =
SRCS     = $(SRCS_F) $(SRCS1_F) $(SRCS_F77) $(SRCS_C)

#
# TARGETs
# all - build the library and install it
all: $(LIB)

$(LIB):$(OBJS) 
	$(AR) $(ARFLAGS) $(EXTRA_ARFLAGS) $@ $?
	$(CP) $(LIB) $(LIBDIR)

clean :
	$(RM) $(LIBDIR)/$(LIB)
	$(RM) $(OBJS) $(LIB)
	@for f in `ls -1 *.F|sed "s/.F$$/.f/"` ; do \
		 $(RM) $$f   ; \
   done

#
# Force serial compiles / preprocessor aligns in free format
.f.o:
	$(SFC) $(FFLAGS) $(EXTRA_FFLAGS) $<

.F.o:
	$(CPP) $(CPP_FLAGS) $< > $*.f
	$(SFC) $(FFLAGS) $(EXTRA_FFLAGS) $*.f

.IGNORE:
.PHONY: clean
