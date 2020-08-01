#
# for the XML-Fortran modules and parsing library
SHELL = /bin/sh -x
LIB   = libxmlparse.a
PGM   = xmlreader

#
#configuration file contains architecture and compile information
include ../../../configure.upp

#
# extra flags
EXTRA_CPPFLAGS  =
EXTRA_FFLAGS    = -c -O2
EXTRA_CFLAGS    = -c
EXTRA_ARFLAGS   =

all             : $(LIB)

#
# PGM xmlreader is a program which creates the template file for
# post to build with - not created or used for building community post
$(PGM)          : xmlreader.o $(LIB)
	$(FC) -o $(PGM) xmlreader.o $(LIB)

$(LIB)          : xmlparse.o read_xml_prims.o
	$(AR) $(ARFLAGS) $@ $?
	$(CP) *.mod  $(INCMOD)
	$(CP) $(LIB) $(LIBDIR) 

clean:
	@for f in `ls *.mod` ; do \
		$(RM) $$f ; $(RM) $(INCMOD)/$$f ; \
	done
	$(RM) *.o
	$(RM) *.mod
	$(RM) *.a
	$(RM) $(LIBDIR)/$(LIB)
	$(RM) $(PGM)

.f90.o:
	$(FC) $(FFLAGS) $(EXTRA_FFLAGS) $<

#-------------Dependencies------------------------------------
read_xml_prims.o : read_xml_array.inc read_xml_scalar.inc read_from_buffer.inc
