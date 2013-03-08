#
# This is a makefile for the UPP libraries
SHELL=/bin/sh 

#
# The configuration file created based on the architecture and compiler
include ../../configure.upp

#
# Note crtm2 library not part of the UPP repository at NCAR - but is available
#   at NCAR in the external-lib repository -- the library is bundled as part
#   of UPP for any official release
SUBDIRS = crtm2 ip sfcio sigio sp w3 nemsio gfsio g2 g2tmpl xml $(SERIAL_MPI_STUB)
LIBDIR  = ../../lib
INCMOD  = ../../include

#----------------------------------------------------------------------------
# bacio must be first other dependencies

all: bacio $(SUBDIRS)
	$(CP) crtm2.makefile crtm2/makefile
	@for dir in bacio $(SUBDIRS); do \
		( cd $$dir; \
		  echo -e "\n===== Making $@ in `pwd` =====" ; \
		  $(MAKE) $@ ) ; \
	done

#-----------------------------------------------------------------------------
# clean all directories -- being tidy

clean: bacio $(SUBDIRS) wrfmpi_stubs
	@for dir in bacio $(SUBDIRS) wrfmpi_stubs; do \
		( cd $$dir; \
		echo -e "\n====== $@ in `pwd` ======" ; \
		$(MAKE) $@ ) ; \
	done
	$(RM) -f crtm2/makefile

.IGNORE:
.PHONY: clean
