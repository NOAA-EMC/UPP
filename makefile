SHELL=/bin/sh

#
# TARGETs

all: $(SUBDIRS)
	cd sorc/comlibs; $(MAKE)
	cd sorc/ncep_post.fd; $(MAKE) -f makefile_dtc

clean: $(SUBDIRS)
	cd sorc/comlibs; $(MAKE) clean
	cd sorc/ncep_post.fd; $(MAKE) clean -f makefile_dtc

.IGNORE:
.PHONY: clean
