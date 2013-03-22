# This makefile compiles copygb
# NOTE: This is a serial build as the code is written to be serial
# Usage: copygb [-g "grid [kgds]"] [-i "ip [ipopts]"] [-k "kpds"] [-v "uparms"]
#              [-B mapgrib [-b mapindex] [-A "<> mapthreshold"] [-K "mapkpds"]]
#              [-M "mask"/mergegrib [-m mergeindex]] [-X] [-a] [-N namelist]
#       then either:
#              grib1 index1 grib2
#            or:
#              -x grib1 grib2
#
# Perform destaggering (NMM) and HORIZONTAL interpolation to a defined output grid.
# Create output grid NOT fixed by model integration domain

SHELL    = /bin/sh
TARGET   = copygb.exe

#
# build configuration based on architecture and compiler
include ../../configure.upp

#
# directories of shared resources
LOCALINC = -I$(INCMOD)

LLIBDIR   = -L$(LIBDIR)
UPPLIBS  = -lip -lbacio -lsp -lw3nco -lw3emc
LIBS     = $(LLIBDIR) $(UPPLIBS)

MODULES  = ../unipost/kinds_mod.o

#
# extra compilation/link flags
#   we are using single word (4 byte) for all libraries
EXTRA_FFLAGS   =
#EXTRA_FFLAGS   = -c -r8
EXTRA_LDFLAGS  = $(LIBS) 

SRCS=	copygb.f
OBJS=   $(SRCS:.f=.o)

#
# TARGETS

all: $(BINDIR)/$(TARGET)

$(BINDIR)/$(TARGET): $(TARGET)
	$(CP) $(TARGET) $(BINDIR)

$(TARGET): $(SRCS) $(MODULES)
	$(SFC) -o $@ $(SRCS) $(MODULES) $(FFLAGS) $(EXTRA_FFLAGS) $(LDFLAGS) $(EXTRA_LDFLAGS)

clean:	
	@ echo -e "\n <><><><> CLEAN <><><><>\n"
	$(RM) $(TARGET) $(OBJS)
	$(RM) $(BINDIR)/$(TARGET)

.IGNORE:
.PHONY: clean
