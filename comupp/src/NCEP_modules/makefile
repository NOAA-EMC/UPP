# NCEP_modules
#
# These are modules that NCEP uses from a shared resource that is not
# available at NCAR - so we add them here
SHELL=/bin/sh

#
# Our configuration that we created before compilation
include ../../configure.upp

#
# Extra flags
EXTRA_FFLAGS = -c $(PRECISION)

#
# Objects to be made
OBJS_F77 =
OBJS_F   = kinds_mod.o constants_mod.o 

OBJS     = $(OBJS_F77) $(OBJS_F)

SRC_F77 =
SRC_F   = $(OBJS:.o=.f)

#
# These are F90 files w/ a .f extension
all:$(OBJS)


# Remove all compilation files 
clean: 
	@echo -e "\n<><><><> CLEAN <><><><>\n$@ in `pwd`"
	@for f in `ls *.mod`; do \
		( $(RM) $$f ; $(RM) $(INCMOD)/$$f ) ; \
	done
	$(RM) $(OBJS)

.IGNORE:
.PHONY: clean

# override for non compliant filename
.f.o:
	$(F90) $(FFLAGS) $(EXTRA_FFLAGS) $(SRC_F)
	$(CP) *.mod $(INCMOD)
