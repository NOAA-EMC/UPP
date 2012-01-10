#-------------------------------------------------------------------------
# This makefile is for the mersenne library
#
# The makefile expects two variables set from the command line
#-------------------------------------------------------------------------

SHELL = /bin/sh
LIB   = libmersenne.a

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
OBJS_F   =
OBJS_F77 = mersenne_twister.o
OBJS_C   =
OBJS     = $(OBJS_F) $(OBJS_F77) $(OBJS_C)

SRCS_F   =
SRCS_F77 = $(OBJS_F77:.0=.f)
SRCS_C   =
SRCS     = $(SRCS_F) $(SRCS_F77) $(SRCS_C)

#
# TARGETs
# all - build the library and install it
all: $(LIB)

$(LIB): $(OBJS)
	$(AR) $(ARFLAGS) $(EXTRA_ARFLAGS) $@ $(OBJS_F77)
	$(CP) $(LIB) $(LIBDIR)
	$(CP) *.mod $(INCMOD)

clean:
	@for f in `ls *.mod` ; do \
		$(RM) $$f ; $(RM) $(INCMOD)/$$f ; \
	done
	$(RM) $(LIBDIR)/$(LIB)
	$(RM) $(OBJS) $(LIB)

.IGNORE:
.PHONY: clean
