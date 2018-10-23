#--------------------------------------------------------------------
# This makefile is for the gfsio library
#
# The makefile expect two variables set from the command line 
#--------------------------------------------------------------------

SHELL = /bin/sh -x
LIB   = libgfsio.a

#
# configuration file contains architecture and compile information
include ../../../configure.upp

#
# extra flags
EXTRA_FFLAGS  = -c $(PROMOTION) -I$(INCMOD)
EXTRA_CFLAGS  = -c 
EXTRA_ARFLAGS =

#
# What files are we looking for - 
OBJS_F   = gfsio_module.o
OBJS_F77 = 
OBJS_C   = 
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
	for f in `ls -1 *.F|sed "s/.F$$/.f/"` ; do \
		$(RM) $$f   ; \
	done
#
# always serial build
.f.o:
	$(SF90) $(FFLAGS) $(EXTRA_FFLAGS) $<
.F.o:
	$(CPP) $(CPP_FLAGS) $< > $*.f
	$(SF90) $(FFLAGS) $(EXTRA_FFLAGS) $*.f


.PHONY: clean
.IGNORE:
