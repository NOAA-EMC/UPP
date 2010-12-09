#
# This is a makefile for the UPP libraries
SHELL=/bin/sh

#
# The configuration file created based on the architecture and compiler
include ../../configure.upp

#
# Note crtm library is a subset of the original retrieved from NCEP
SUBDIRS = upp_crtm ip mersenne sfcio sigio sp w3 
LIBDIR  = ../../lib
INCMOD  = ../../include

#----------------------------------------------------------------------------
# bacio must be first other dependencies

all: bacio $(SUBDIRS)
	@for dir in bacio $(SUBDIRS); do \
		( cd $$dir; \
		  echo -e "\n===== Making $@ in `pwd` =====" ; \
		  make ) ; \
	done

clean: bacio $(SUBDIRS)
	@for dir in bacio $(SUBDIRS); do \
      ( cd $$dir; \
        echo -e "\n====== $@ in `pwd` ======" ; \
        make $@ ) ; \
   done
