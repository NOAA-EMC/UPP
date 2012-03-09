#
# for the XML-Fortran module
#
FC      =       ncepxlf90
FOPT    =       -c -O2

all             :       libxmlparse.a xmlreader

install         :
	cp libxmlparse.a ../../../xmlfortran/.
	cp xmlparse.mod ../../../xmlfortran/incmod

libxmlparse.a   :       xmlparse.o read_xml_prims.o
	ar r libxmlparse.a xmlparse.o read_xml_prims.o

xmlparse.o      :       xmlparse.f90
	$(FC) $(FOPT) xmlparse.f90

xmlreader       :       xmlreader.f90 libxmlparse.a
	$(FC) $(FOPT) xmlreader.f90
	$(FC) -o xmlreader xmlreader.o libxmlparse.a

read_xml_prims.o        :       read_xml_prims.f90 xmlparse.f90
	$(FC) $(FOPT) read_xml_prims.f90

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.a
