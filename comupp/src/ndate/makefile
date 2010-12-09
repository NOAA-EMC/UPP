# ndate utility program 
# Usage: ndate [fhour [idate]]
#  give the date in the format yyyymmddhh based on the current
#  time w/o arguments or calc to output format a forcast hour 
#  and initial time
SHELL      = /bin/sh
TARGET     = ndate

#
# Configuration created based on architecture and compiler
include ../../configure.upp

#
# Install directory for executable
BINDIR     = ../../bin

#
# Shared resources
LIBDIR     = ../../lib
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
	$(FC) -o $@ $(LDFLAGS) $(EXTRA_LDFLAGS) $(LIBS) $(OBJ)
	$(CP) $(TARGET) ../../bin

clean:
	@ echo -e "\n <><><><> CLEAN <><><><>\n"
	@$(RM) $(TARGET) $(OBJ)
	@$(RM) $(BINDIR)/$(TARGET)
