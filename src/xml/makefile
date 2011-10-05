#LIB     = ./libxmlread.a
INCMOD  = ./incmod/xmlread
FC      =       ncepxlf90
FOPT    =       -g -c -I ../Grumbine/xml-fortran/src
AR      = ar
ARFLAGS = -rv -X64

all: libxmlread.a

libxmlread.a        :   postcntrl_t.o
	ar -ruv libxmlread.a  postcntrl_t.o
	mv xml_data_postcntrl_t.mod $(INCMOD)

#postcntrl_t.o      :       postcntrl_t.f90
#	$(FC) $(FOPT) postcntrl_t.f90

postcntrl_t.f90    :       postcntrl_t.xml
	echo postcntrl_t > xmlreader.inp
	../Grumbine/xml-fortran/src/xmlreader
	$(FC) $(FOPT) postcntrl_t.f90

clean:
	 rm -f postcntrl_t.f90 *.o
	 rm -f *.mod
	 rm -f *.a
