#-------------------------------------------------------------------------
# This makefile is for the sigio library
#
# The makefile expects two variables set from the command line
#-------------------------------------------------------------------------

SHELL = /bin/sh
LIB   = libsigio.a

#
# configuration based on architecture and compiler
include ../../../configure.upp

#
#extra flags
EXTRA_FFLAGS  = -c $(PROMOTION) -I$(INCMOD)
EXTRA_CFLAGS  =
EXTRA_ARFLAGS =

#
# Library object files
OBJS_F   = sigio_module.o sigio_r_module.o
OBJS_F77 = 
OBJS_C   =
OBJS     = $(OBJS_F) $(OBJS_F77) $(OBJS_C)

#
# Source files - non conventional naming
SRCS_F   = $(OBJS_F:.o=.f)
SRCS_F77 = $(OBJS_F77:.o=.f)
SRCS_C   = $(OBJS_C:.o=.c)
SRCS     = $(SRCS_F) $(SRCS_F77) $(SRCS_C)

all: $(LIB)

$(LIB): $(OBJS_F)
	$(AR) $(ARFLAGS) $(EXTRA_ARFLAGS) $@ $?
	$(CP) $(LIB) $(LIBDIR)
	$(CP) *.mod $(INCMOD)

clean:
	@for f in `ls *.mod` ; do \
		$(RM) $$f ; $(RM) $(INCMOD)/$$f ; \
	done
	$(RM) $(LIBDIR)/$(LIB)
	$(RM) $(OBJS_F) *.mod $(LIB)
	$(RM) ${INCMOD}/sigio_module.mod ${INCMOD}/sigio_r_module.mod

.IGNORE:
.PHONY: clean

# override non convention compliant filename
.SUFFIXES : .f
.f.o:
	$(F90) $(FFLAGS) $(EXTRA_FFLAGS) $<
