#
# This is transitional makefile from Community UPP to the CRTM library
#
#==============================================================================
#
# Community Radiate Transfer Model (CRTM) Makefile
#   The library name is set in the CRTM makefiles -- so this file must
#   be modified if CRTM uodates theirs
#==============================================================================

SHELL = /bin/sh
LIB   = libCRTM.a
INCMOD_CRTM = $(INCMOD)/crtm2

#
# configuration file contains architecture and compile information
include ../../../configure.upp

#
# Needed system commands

#
# Extra Flags
EXTRA_FFLAGS  = -c $(PROMOTION)
EXTRA_CFLAGS  = -c
EXTRA_ARFLAGS =

#
# TARGETS
all :
	( cd src; echo "Making CRTM library in `pwd`" ; \
	  $(MAKE) FC="$(F90)" FL="$(F90)" FC_FLAGS="$(FFLAGS_CRTM)" FL_FLAGS="$(FL_CRTM)" install; \
		\
	  $(CP) lib/libCRTM.a $(LIBDIR)/$(LIB) ; \
	  $(LN) `pwd`/include $(INCMOD_CRTM) ; \
	)

#
# Make clean - always use crtm distclean
clean:
	( cd src ; echo "Cleaning CRTM library" ; \
	  $(MAKE) distclean ; \
	  $(RM) $(INCMOD_CRTM) ; \
	  $(RM) $(LIBDIR)/$(LIB) ; \
   )
distclean: clean

.IGNORE:
.PHONY: distclean clean
