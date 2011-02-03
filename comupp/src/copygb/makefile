# This makefile compiles copygb
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

SHELL    =	/bin/sh
TARGET   = copygb.exe

#
# build configuration based on architecture and compiler
include ../../configure.upp

#
# directories of shared resources
LOCALINC = -I$(INCMOD)

LLIBDIR   = -L$(LIBDIR)
UPPLIBS  = -lw3 -lip -lbacio -lsp	
LIBS     = $(LLIBDIR) $(UPPLIBS)

MODULES  = ../NCEP_modules/kinds_mod.o ../NCEP_modules/constants_mod.o

#
# extra compilation/link flags
EXTRA_FFLAGS   = -c
#EXTRA_FFLAGS   = -c -r8
EXTRA_LDFLAGS  = $(LIBS) 


OBJS=	copygb.o

#
# TARGETS

################################################################################
all: $(TARGET)

$(TARGET):	$(OBJS)
	$(FC) -o $@ $(OBJS) $(MODULES) $(LDFLAGS) $(EXTRA_LDFLAGS)
	$(CP) $(TARGET) $(BINDIR)

clean:	
	@ echo -e "\n <><><><> CLEAN <><><><>\n"
	@$(RM) $(TARGET) $(OBJS)
	@$(RM) $(BINDIR)/$(TARGET)
