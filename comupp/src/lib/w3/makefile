#--------------------------------------------------------------------
# This makefile is for the w3 library
#
# The makefile expect two variables set from the command line
#--------------------------------------------------------------------

SHELL = /bin/sh
LIB   = libw3.a

#
# build configuration determined before compile
include ../../../configure.upp

#
# Extra Flags
EXTRA_FFLAGS  = -c $(PROMOTION) -I$(INCMOD)
EXTRA_CFLAGS  = -c
EXTRA_ARFLAGS =

#
# What files are we expecting - note there is a FORTRAN mova2i and C -- what's
#   up with that?
OBJS_F   =
OBJS_F77	= w3fa09.o w3fb06.o w3fb11.o errexit.o errmsg.o fparsei.o fparser.o \
           gbytes_char.o getbit.o getgb1re.o getgb1r.o getgb1s.o    \
           getgbe.o getgbeh.o getgbem.o getgbemh.o getgbemp.o getgbens.o     \
           getgbep.o getgbex.o getgbexm.o getgb.o getgbh.o getgbm.o getgbmh.o\
           getgbmp.o getgbp.o getgi.o getgir.o idsdef.o instrument.o iw3jdn.o\
           ixgb.o lengds.o mova2i.o pdsens.o pdseup.o putgbe.o putgben.o     \
           putgbens.o putgbex.o putgb.o putgbn.o r63w72.o    \
           skgb.o start.o summary.o w3valdat.o w3difdat.o w3doxdat.o w3fi01.o\
           w3fi58.o w3fi59.o w3fi63.o w3fi68.o w3fi71.o w3fi72.o w3fi73.o    \
           w3fi74.o w3fi75.o w3fi76.o w3fi82.o w3fi83.o w3fs21.o w3fs26.o    \
           w3locdat.o w3log.o w3movdat.o w3reddat.o w3tagb.o w3utcdat.o      \
           mersenne_twister.o
OBJS_C   = mova2i.o

OBJS     = $(OBJS_F) $(OBJS_F77) $(OBJS_C)

#
# Source files
SRCS_F   = $(OBJS_F:.o=.f90)
SRCS_F77 = $(OBJS_F77:.o=.f)
SRCS_C   = $(OBJS_C:.o=.c)
SRCS     = $(SRCS_F) $(SRCS_F77) $(SRCS_C)

#
# TARGETs
# all - build the library and install it
all: $(LIB)

$(LIB): $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(CP) $(LIB) $(LIBDIR)
	$(CP) *.mod $(INCMOD)

clean:
	@for f in `ls *.mod` ; do \
		$(RM) $$f ; $(RM) $(INCMOD)/$$f ; \
	done
	$(RM) $(LIBDIR)/$(LIB)
	$(RM) $(OBJS) $(LIB) *.mod

# Always serial build
.f.o:
	$(SFC) $(FFLAGS) $(EXTRA_FFLAGS) $<

.IGNORE:
.PHONY: clean
