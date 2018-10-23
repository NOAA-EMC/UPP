SHELL=/bin/sh
LIB=libg2.a

#
# configuration file contains architecture and compile information
include ../../../configure.upp

#  If you want to enable support for PNG or JPEG2000 encoding/decoding,
#  you must specify -DUSE_PNG and/or -DUSE_JPEG2000 in the FDEFS variable
#  for the Fortran pre-processor
#  -DUSE_PNG requires libpng.a and libz.a
#                     ( and png.h pngconf.h zconf.h zlib.h include files).
#  -DUSE_JPEG2000 requires libjasper.a
#                     ( and all the jasper/*.h include files).
#
#  In addition, INCDIR must include all directories where the above
#  mentioned include files can be found.
FDEFS=-DUSE_JPEG2000

#
# extra flags
EXTRA_CPPFLAGS = $(FDEFS)
EXTRA_FFLAGS  = -c $(PROMOTION) -I$(INCMOD)
EXTRA_CFLAGS  = -c $(GRIB2SUPT_INC)
EXTRA_ARFLAGS =
 
#--------------------------------------
# TARGETs

all: $(LIB)
	$(CP) $(LIB) $(LIBDIR)
	$(CP) *.mod $(INCMOD)

$(LIB):	$(LIB)(gridtemplates.o) \
	$(LIB)(pdstemplates.o) \
	$(LIB)(drstemplates.o) \
	$(LIB)(gribmod.o) \
	$(LIB)(realloc.o) \
	$(LIB)(addfield.o) \
	$(LIB)(addgrid.o) \
	$(LIB)(addlocal.o) \
	$(LIB)(getfield.o) \
	$(LIB)(gb_info.o) \
	$(LIB)(gf_getfld.o) \
	$(LIB)(gf_free.o) \
	$(LIB)(gf_unpack1.o) \
	$(LIB)(gf_unpack2.o) \
	$(LIB)(gf_unpack3.o) \
	$(LIB)(gf_unpack4.o) \
	$(LIB)(gf_unpack5.o) \
	$(LIB)(gf_unpack6.o) \
	$(LIB)(gf_unpack7.o) \
	$(LIB)(gettemplates.o) \
	$(LIB)(getlocal.o) \
	$(LIB)(getdim.o) \
	$(LIB)(getpoly.o) \
	$(LIB)(gribcreate.o) \
	$(LIB)(gribend.o) \
	$(LIB)(gribinfo.o) \
	$(LIB)(mkieee.o) \
	$(LIB)(rdieee.o) \
	$(LIB)(simpack.o) \
	$(LIB)(simunpack.o) \
	$(LIB)(cmplxpack.o) \
	$(LIB)(compack.o) \
	$(LIB)(misspack.o) \
	$(LIB)(pack_gp.o) \
	$(LIB)(reduce.o) \
	$(LIB)(comunpack.o) \
	$(LIB)(specpack.o) \
	$(LIB)(specunpack.o) \
	$(LIB)(jpcpack.o) \
	$(LIB)(jpcunpack.o) \
	$(LIB)(enc_jpeg2000.o) \
	$(LIB)(dec_jpeg2000.o) \
	$(LIB)(pngpack.o) \
	$(LIB)(pngunpack.o) \
	$(LIB)(enc_png.o) \
	$(LIB)(dec_png.o) \
        $(LIB)(mova2i.o) \
	$(LIB)(gbytesc.o) \
	$(LIB)(skgb.o) \
	$(LIB)(ixgb2.o) \
        $(LIB)(getidx.o) \
	$(LIB)(getidxes.o)\
	$(LIB)(getg2i.o) \
	$(LIB)(getg2ir.o) \
	$(LIB)(getgb2s.o) \
	$(LIB)(getgb2r.o) \
	$(LIB)(getgb2l.o) \
	$(LIB)(getgb2.o) \
	$(LIB)(getgb2p.o) \
	$(LIB)(getgb2rp.o) \
	$(LIB)(putgb2.o) \
	$(LIB)(g2grids.o) \
        $(LIB)(gdt2gds.o)\
	$(LIB)(params.o) \
        $(LIB)(intmath.o)\
        $(LIB)(params_ecmwf.o)

$(LIB)(compack.o): $(LIB)(intmath.o)
$(LIB)(simpack.o): $(LIB)(intmath.o)
$(LIB)(misspack.o): $(LIB)(intmath.o)

clean :
	@for f in `ls *.mod` ; do \
		$(RM) $$f ; $(RM) $(INCMOD)/$$f ; \
	done
	@for f in `ls -1 *.F|sed "s/.F$$/.f/"` ; do \
		$(RM) $$f   ; \
	done
	$(RM) $(LIBDIR)/$(LIB)
	$(RM) $(LIB)
	$(RM) $*.o

.F.a:
	$(CPP) $(CPP_FLAGS) $(EXTRA_CPPFLAGS) $*.F $*.f
	$(FC) -c $(FFLAGS) $(EXTRA_FFLAGS) $*.f
	$(AR) $(ARFLAGS) $(EXTRA_ARFLAGS) $@ $*.o

.f.a:
	$(FC) -c $(FFLAGS) $(EXTRA_FFLAGS) $<
	$(AR) $(ARFLAGS) $(EXTRA_ARFLAGS) $@ $*.o

.c.a:
	$(CC) -c $(CFLAGS) $(EXTRA_CFLAGS) $<
	$(AR) $(ARFLAGS) $(EXTRA_ARFLAGS) $@ $*.o
