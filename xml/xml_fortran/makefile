# Naive makefile for the simple test program "tstparse"
#
# $Id: makefile,v 1.6 2006/06/09 10:30:30 arjenmarkus Exp $
#
# Coarse test for the XML-Fortran module
#
# Note:
# You may need to change the name of the Fortran compiler
# and the compile options
#
#FC	=	f90
#FC	=	lf95
#FC	=	ifort
#FOPT	=	-c
FC	=	ncepxlf90
FOPT	=	-c -O2

all		:	xmlparse.a tstparse xmlreader

install		:	
	cp xmlparse.a ../lib
	cp xmlparse.mod ../lib

xmlparse.a	:	xmlparse.o read_xml_prims.o
	ar r xmlparse.a xmlparse.o read_xml_prims.o

xmlparse.o	:	xmlparse.f90
	$(FC) $(FOPT) xmlparse.f90

tstparse	:	tstparse.f90 xmlparse.a
	$(FC) $(FOPT) tstparse.f90
	$(FC) -o tstparse tstparse.o xmlparse.a

xmlreader	:	xmlreader.f90 xmlparse.a
	$(FC) $(FOPT) xmlreader.f90
	$(FC) -o xmlreader xmlreader.o xmlparse.a

read_xml_prims.o	:	read_xml_prims.f90 xmlparse.f90
	$(FC) $(FOPT) read_xml_prims.f90

clean:	
	rm -f *.o
	rm -f *.mod
	rm -f *.a
