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
LIBS       = $(LIBDIR)/libw3.a

#
# Extra flags
EXTRA_FFLAGS  = -c

#
# Objects to be made
OBJ        =   ndate.o

#
# TARGETS
all: $(TARGET)

$(TARGET): $(OBJ)
	$(SFC) -o $@ $(OBJ) $(LDFLAGS) $(EXTRA_LDFLAGS) $(LIBS)
	$(CP) $(TARGET) ../../bin

clean:
	@ echo -e "\n <><><><> CLEAN <><><><>\n"
	$(RM) $(TARGET) $(OBJ)
	$(RM) $(BINDIR)/$(TARGET)

.IGNORE:
.PHONY: clean
