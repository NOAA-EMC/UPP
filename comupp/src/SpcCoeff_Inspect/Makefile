#==============================================================================
#
# Makefile for SpcCoeff_Inspect program
#
#==============================================================================

# Executable files
EXE_FILE = ../../bin/spccoeff_inspect.exe

all: $(EXE_FILE)

clean: 
	rm -f $(EXE_FILE) *.o a.out *~ *.mod

# Define macros
CRTM_SOURCE_ROOT=../lib/crtm2/src


include  ../../configure.upp
include $(CRTM_SOURCE_ROOT)/make.macros

# This makefile
MAKE_FILE = Makefile


$(EXE_FILE): SpcCoeff_Inspect.f90
	$(SF90) -o $(EXE_FILE) $(FFLAGS) -I$(CRTM_SOURCE_ROOT)/libsrc SpcCoeff_Inspect.f90 -L../../lib/ -lCRTM

# Source files to link
SRC_FILES = Type_Kinds.f90 \
            File_Utility.f90 \
            Message_Handler.f90 \
            Compare_Float_Numbers.f90 \
            Binary_File_Utility.f90 \
            Endian_Utility.f90 \
            Sort_Utility.f90 \
            Subset_Define.f90 \
            SensorInfo_Parameters.f90 \
            SpcCoeff_Define.f90 \
              ACCoeff_Define.f90 \
              NLTECoeff_Define.f90 \
            SpcCoeff_Binary_IO.f90 \
              ACCoeff_Binary_IO.f90 \
              NLTECoeff_Binary_IO.f90

# Obj files used in link phase
OBJ_FILES = ${SRC_FILES:.f90=.o} \
            $(EXE_FILE).o

# Define common make targets (all, build, clean, install)
#include $(CRTM_SOURCE_ROOT)/make.common_targets

# Source link creation and removal
#create_links::
#	@$(LINK_SCRIPT) $(CRTM_SOURCE_ROOT) $(SRC_FILES)

#remove_links::
#	@$(UNLINK_SCRIPT) $(SRC_FILES)

# Source dependency lists
#include make.dependencies

# Define default rules
include $(CRTM_SOURCE_ROOT)/make.rules

