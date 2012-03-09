LIB     = ./libxmlread.a
INCMOD  = ./xml_fortan
FC      =       ncepxlf90
FOPT    =       -g -c -I ./xml_fortran
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
	./xml_fortran/xmlreader
	$(FC) $(FOPT) postcntrl_t.f90

clean:
	 rm -f postcntrl_t.f90 *.o
	 rm -f *.mod
	 rm -f *.a
