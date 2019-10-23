#--------------------------------------------------------------------
# This makefile is for the w3 library
#
# The makefile expect two variables set from the command line
#--------------------------------------------------------------------

SHELL = /bin/sh
LIB   = libw3emc.a

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
OBJS_F77 = w3fa09.o w3fb06.o w3fb11.o getgbens.o     \
           mova2i.o putgbens.o start.o summary.o     \
           w3log.o mersenne_twister.o
OBJS_C   = 

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
