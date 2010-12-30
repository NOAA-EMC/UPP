#-------------------------------------------------------------------------
# This makefile is for the sfcio library
#
# The makefile expects two variables set from the command line
#-------------------------------------------------------------------------

SHELL = /bin/sh
LIB   = libsfcio.a

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
OBJS_F = sfcio_module.o

SRCS_F = $(OBJS_F:.0=.f)

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
	$(RM) $(OBJS_F) $(LIB)

# Override non convention compliant filename
.SUFFIXES : .f
.f.o:
	$(F90) $(FFLAGS) $(EXTRA_FFLAGS) $<
