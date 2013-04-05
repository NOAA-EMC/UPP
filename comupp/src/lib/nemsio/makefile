#--------------------------------------------------------------------
# This makefile is for the bacio library
#
# The makefile expect two variables set from the command line 
#--------------------------------------------------------------------

SHELL = /bin/sh -x
LIB   = libnemsio.a

#
# configuration file contains architecture and compile information
include ../../../configure.upp

#
# extra flags
EXTRA_FFLAGS  = -c $(PROMOTION) -I$(INCMOD)
EXTRA_CFLAGS  = -c 
EXTRA_ARFLAGS =

#
# What files are we looking for - yes bafrio, very poor style
OBJS_F   = nemsio_openclose.o nemsio_read.o nemsio_write.o \
           nemsio_module.o bafrio.o
OBJS_F77 =
OBJS_C   = 
OBJS     = $(OBJS_F) $(OBJS_F77) $(OBJS_C)

SRCS_F   = $(OBJS_F:.o=.f90)
SRCS_F77 = $(OBJS_F77:.o=.f)
SRCS_C   = $(OBJS_C:.o=.c)
SRCS     = $(SRCS_F) $(SRCS_F77) $(SRCS_C)

#
# TARGETs
# all - build the library and install it
all: $(LIB)

$(LIB): $(OBJS) $(OBJS_F77) $(OBJS_C)
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
.f90.o:
	$(F90) $(FFLAGS) $(EXTRA_FFLAGS) $<

#
# file dependencies
nemsio_read.o: nemsio_openclose.o
nemsio_write.o: nemsio_openclose.o
nemio_module.o: nemsio_openclose.o nemsio_read.o nemsio_write.o

.PHONY: clean
.IGNORE:
