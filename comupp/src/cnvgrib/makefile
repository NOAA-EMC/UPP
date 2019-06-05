# ndate utility program 
# Usage: ndate [fhour [idate]]
#  give the date in the format yyyymmddhh based on the current
#  time w/o arguments or calc to output format a forcast hour 
#  and initial time
SHELL      = /bin/sh
TARGET     = cnvgrib.exe

#
# Configuration created based on architecture and compiler
include ../../configure.upp

#libbacio.a  libCRTM.a  libg2.a  libg2tmpl.a  libgfsio.a  libip.a  libnemsio.a  libsfcio.a  libsigio.a  libsp.a  libw3emc.a  libw3nco.a  libxmlparse.a

#
# Shared resources
LIBS       = -L$(LIBDIR) -lg2 -lw3nco -lbacio $(GRIB2SUPT_LIB)
FFLAGS += -I../lib/g2  # for g2/*.mod files

#
# Extra flags
EXTRA_FFLAGS  =

#
# Objects to be made
SRC=	gds2gdt.f pds2pdt.f pds2pdtens.f cnvgrib.f cnv12.f cnv22.f \
        cnv21.f gdt2gds.f makepds.f putgbexn.f makepdsens.f setbit.f

OBJ        = $(SRC:.f=.o)

#
# TARGETS
all: $(BINDIR)/$(TARGET)

$(BINDIR)/$(TARGET): $(TARGET)
	$(CP) $(TARGET) $(BINDIR)

$(TARGET): $(SRC)
	$(SFC) -o $@ $(SRC) $(FFLAGS) $(EXTRA_FFLAGS) $(LDFLAGS) $(EXTRA_LDFLAGS) $(LIBS)

clean:
	@ echo -e "\n <><><><> CLEAN <><><><>\n"
	$(RM) $(TARGET) $(OBJ)
	$(RM) $(BINDIR)/$(TARGET)

.IGNORE:
.PHONY: clean
