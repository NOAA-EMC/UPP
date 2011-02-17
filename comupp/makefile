SHELL=/bin/sh

#
# Bring in configuration determined for UPP
#include ./configure.upp

#
# dependent build order of source -
# unipost relies on NCEP_modules and lib builds
# copygb  relies on NCEP_modules and lib builds
# ndate   relies on lib builds
SUBDIRS = src/NCEP_modules src/lib src/unipost src/copygb src/ndate

#
# TARGETs

all: $(SUBDIRS)
	@for dir in $(SUBDIRS); do \
      ( cd $$dir; echo "Making $@ in `pwd`" ; \
        make $@ ); \
   done

clean: $(SUBDIRS)
	@for dir in $(SUBDIRS); do \
      ( cd $$dir; echo "Making $@ in `pwd`" ; \
        make $@) ; \
   done

.IGNORE
.PHONY: clean
