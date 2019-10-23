#--------------------------------------------------------------------
# This makefile is for the bacio library
#
# The makefile expect two variables set from the command line 
#--------------------------------------------------------------------

SHELL = /bin/sh -x
LIB   = libbacio.a

#
# configuration file contains architecture and compile information
include ../../../../configure.upp

#
# extra flags
EXTRA_FFLAGS  = -c $(PROMOTION) -I$(INCMOD)
EXTRA_CFLAGS  = -c 
EXTRA_ARFLAGS =

#
# What files are we looking for - note the convention for noting the free
# format of f90 was not maintained so we account for it manually (sorry)
OBJS_F   = 
OBJS_F77 = bafrio.o baciof.o chk_endianc.o
OBJS_C   = bacio.o byteswap.o
OBJS     = $(OBJS_F) $(OBJS_F77) $(OBJS_C)

SRCS_F   = $(OBJS_F:.o=.f)
SRCS_F77 = $(OBJS_F77:.o=.f)
SRCS_C   = $(OBJS_C:.o=.c)
SRCS     = $(SRCS_F) $(SRCS_F77) $(SRCS_C)

#
# TARGETs
# all - build the library and install it
all: $(LIB)

$(LIB): $(OBJS) 
	$(AR) $(ARFLAGS) $(EXTRA_ARFLAGS) $@ $(OBJS)
	$(CP) $(LIB) $(LIBDIR)
	$(CP) *.mod $(INCMOD)

clean :
	@for f in `ls *.mod` ; do \
		$(RM) $$f ; $(RM) $(INCMOD)/$$f ; \
	done
	$(RM) $(LIBDIR)/$(LIB)
	$(RM) $(OBJS_F) $(OBJS_F77) $(OBJS_C) $(LIB)
#
# always serial build
.f.o:
	$(SFC) $(FFLAGS) $(EXTRA_FFLAGS) $<

.PHONY: clean
.IGNORE:
