#--------------------------------------------------------------------
# This makefile is for the bacio library
#
# The makefile expect two variables set from the command line 
#--------------------------------------------------------------------

SHELL = /bin/sh
LIB   = libbacio.a

#
# configuration file contains architecture and compile information
include ../../../configure.upp

#
# extra flags
EXTRA_FFLAGS  = -c $(PROMOTION) -I$(INCMOD)
EXTRA_CFLAGS  = -c -Dfunder -DFortranByte=char -DFortranInt=int -DFortranLlong='long long'
EXTRA_ARFLAGS =

#
# What files are we looking for - note the convention for noting the free
# format of f90 was not maintained so we account for it manually (sorry)
OBJS_F   = bafrio.o
OBJS_F77 = baciof.o
OBJS_C   = bacio.v1.4.o
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
	$(AR) $(ARFLAGS) $(EXTRA_ARFLAGS) $@ $?
	$(CP) $(LIB) $(LIBDIR)
	$(CP) *.mod $(INCMOD)

#
# override for non compliant files
$(OBJS_F): $(SRCS_F)
	$(F90) $(FFLAGS) $(EXTRA_FFLAGS)  $<

clean :
	@for f in `ls *.mod` ; do \
		$(RM) $$f ; $(RM) $(INCMOD)/$$f ; \
	done
	$(RM) $(LIBDIR)/$(LIB)
	$(RM) $(OBJS_F) $(OBJS_F77) $(OBJS_C) $(LIB)
