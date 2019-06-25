
include ../../configure.upp

# Thse four lines are copied straight from ../unipost/Makefile:
LLIBDIR     = -L$(LIBDIR)
GRIB2LIBS   = -lg2 -lg2tmpl -lxmlparse $(GRIB2SUPT_LIB)
UPPLIBS     = -lnemsio -lsigio -lsfcio -lgfsio -lsp -lw3nco -lw3emc -lbacio -lCRTM $(SERIAL_MPI_LIB)
LIBS        = $(LLIBDIR) $(UPPLIBS) $(GRIB2LIBS) $(WRFEXTLIBS) $(NCDFLIBS)

EXE=../../bin/satgrib2.exe

all: $(EXE)

clean:
	rm -f *.o *.mod satgrib2 a.out $(EXE) *~

$(EXE): satgrib2.F90 gds2gdt.f
	$(SFC) $(FFLAGS)      \
		 -o $(EXE) satgrib2.F90 gds2gdt.f $(LIBS)
