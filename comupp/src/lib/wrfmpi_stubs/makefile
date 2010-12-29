#-------------------------------------------------------------------------
# This makefile is for the wrfmpi_stubs library
#-------------------------------------------------------------------------

SHELL = /bin/sh
LIB   = libmpi.a

#
# configuration based on architecture and compiler
include ../../../configure.upp

#
# directories for shared resources
INCMOD  = ../../../include
LIBDIR  = ../../../lib

#
#extra flags
EXTRA_FFLAGS  = -c $(PROMOTION) -I$(INCMOD)
EXTRA_CFLAGS  = -c $(PROMOTION) -I$(INCMOD)
EXTRA_ARFLAGS =

#
# Library object files
OBJS_F   =
OBJS_F77 = mpi_fortran.o
OBJS_C   = mpi_c.o
OBJS     = $(OBJS_F) $(OBJS_F77) $(OBJS_C)

SRCS_F   =
SRCS_F77 = $(OBJS_F77:.o=.f)
SRCS_C   = $(OBJS_C:.o=.c)
SRCS     = $(SRCS_F) $(SRCS_F77) $(SRCS_C)

INCLUDE  = mpi.h

#
# TARGETs
# all - build the library and install it
all: $(LIB)

#
# Link must be made when using mpi stub for serial compiles
$(LIB): $(OBJS)
	$(AR) $(ARFLAGS) $(EXTRA_ARFLAGS) $@ $(OBJS)
	$(CP) $(LIB) $(LIBDIR)
	if [ $(SERIAL_MPI_STUB) != "" ] ; then \
		$(LN) ../lib/wrfmpi_stubs/mpif.h ../../unipost/mpif.h ; \
	fi

clean:
	$(RM) ../../unipost/mpif.h
	$(RM) $(LIBDIR)/$(LIB)
	$(RM) $(OBJS) $(LIB)
