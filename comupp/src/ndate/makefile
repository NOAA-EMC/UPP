# ndate utility program 
# Usage: ndate [fhour [idate]]
#  give the date in the format yyyymmddhh based on the current
#  time w/o arguments or calc to output format a forcast hour 
#  and initial time
SHELL      = /bin/sh
TARGET     = ndate.exe

#
# Configuration created based on architecture and compiler
include ../../configure.upp

#
# Shared resources
LIBS       = $(LIBDIR)/libw3nco.a

#
# Extra flags
EXTRA_FFLAGS  =

#
# Objects to be made
SRC        =   ndate.f

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
