#--------------------------------------------------------------------
# This makefile is for the ip library
#
# The makefile expects two variables set from the command line
#--------------------------------------------------------------------

SHELL = /bin/sh
LIB   = libip.a

#
# configuration based on architecture and compiler
include ../../../configure.upp

#
# extra flags
EXTRA_FFLAGS  = -c $(PROMOTION) -I$(INCMOD)
EXTRA_CFLAGS  =
EXTRA_ARFLAGS =

#
# Object files for this library - all are o
OBJS_F   =  gausslat.o \
 gdswzd00_mod.o \
 gdswzd01_mod.o \
 gdswzd03_mod.o \
 gdswzd04_mod.o \
 gdswzd05_mod.o \
 gdswzdcb_mod.o \
 gdswzdcd_mod.o \
 gdswzd_mod.o \
 ijkgds0.o \
 ijkgds1.o \
 ipolates.o \
 ipolatev.o \
 ipxwafs2.o \
 ipxwafs3.o \
 ipxwafs.o \
 makgds.o \
 movect.o \
 polates0.o \
 polates1.o \
 polates2.o \
 polates3.o \
 polates4.o \
 polates6.o \
 polatev0.o \
 polatev1.o \
 polatev2.o \
 polatev3.o \
 polatev4.o \
 polatev6.o \
 polfixs.o \
 polfixv.o

#OBJS_F77 =                                    \
  gausslat.o gcdist.o gdsawt.o gdswiz00.o     \
  gdswiz01.o gdswiz03.o gdswiz04.o gdswiz05.o \
  gdswizc9.o gdswizca.o gdswizcb.o gdswizcd.o \
  gdswiz.o gdswzd00.o gdswzd01.o gdswzd03.o   \
  gdswzd04.o gdswzd05.o gdswzdc9.o gdswzdca.o \
  gdswzdcb.o gdswzdcd.o gdswzd.o ijkgds0.o    \
  ijkgds1.o ijkgds.o ipmerge2.o ipolates.o    \
  ipolatev.o ipsector.o ipspaste.o ipxetas.o  \
  ipxwafs2.o ipxwafs3.o ipxwafs.o makgds.o    \
  movect.o polateg0.o polateg1.o polateg4.o   \
  polates0.o polates1.o polates2.o polates3.o \
  polates4.o polates6.o polatev0.o polatev1.o \
  polatev2.o polatev3.o polatev4.o polatev6.o \
  polfixs.o polfixv.o
OBJS_C   =
OBJS     = $(OBJS_F) $(OBJS_F77) $(OBJS_C)

#
# Source files
SRCS_F   = $(OBJS_F:.o=.f90)
SRCS_F77 = $(OBJS_F77:.o=.f)
SRCS_C   = $(OBJS_C:.o=.c)
SRCS     = $(SRCS_F) $(SRCS_F77) $(SRCS_C)

#
# TARGETS
all: $(LIB)

#---------------------
# Main compilation and archiving of the library

$(LIB): $(OBJS)
	$(AR) $(ARFLAGS) $(EXTRA_ARFLAGS) $@ $?
	$(CP) $(LIB) $(LIBDIR)
#$(CP) *.mod $(INCMOD)

clean:
	@for f in `ls *.mod` ; do \
                $(RM) $$f ; $(RM) $(INCMOD)/$$f ; \
        done
	$(RM) $(LIBDIR)/$(LIB)
	$(RM) $(OBJS) $(LIB)

# Always serial build
.f.o:
	$(SFC) $(FFLAGS) $(EXTRA_FFLAGS) $<
.f90.o:
	$(F90) $(FFLAGS) $(EXTRA_FFLAGS) $<
.F90.o:
	$(CPP) $(CPP_FLAGS) $(EXTRA_CPPFLAGS) $*.F90 $*.f
	$(F90) $(FFLAGS) $(EXTRA_FFLAGS) $<

# file dependencies
gdswzd_c.o : gdswzd_mod.o gdswzd_c.F90
gausslat.o : gausslat.f90
gdswzd00_mod.o : gdswzd00_mod.f90
gdswzd01_mod.o : gdswzd01_mod.f90
gdswzd03_mod.o : gdswzd03_mod.f90
gdswzd04_mod.o : gdswzd04_mod.f90
gdswzd05_mod.o : gdswzd05_mod.f90
gdswzdcb_mod.o : gdswzdcb_mod.f90
gdswzdcd_mod.o : gdswzdcd_mod.f90
gdswzd_mod.o : gdswzd00_mod.o gdswzd01_mod.o gdswzd03_mod.o gdswzd04_mod.o gdswzd05_mod.o gdswzdcb_mod.o gdswzdcd_mod.o gdswzd_mod.f90
ijkgds0.o : ijkgds0.f90
ijkgds1.o : ijkgds1.f90
ipolates.o : ipolates.f90
ipolatev.o : ipolatev.f90
ipxwafs2.o : ipxwafs2.f90
ipxwafs3.o : ipxwafs3.f90
ipxwafs.o : ipxwafs.f90
makgds.o : makgds.f90
movect.o : movect.f90
polates0.o : polates0.f90
polates1.o : polates1.f90
polates2.o : polates2.f90
polates3.o : polates3.f90
polates4.o : polates4.f90
polates6.o : polates6.f90
polatev0.o : polatev0.f90
polatev1.o : polatev1.f90
polatev2.o : polatev2.f90
polatev3.o : polatev3.f90
polatev4.o : polatev4.f90
polatev6.o : polatev6.f90
polfixs.o : polfixs.f90
polfix.o : polfixv.f90

.IGNORE:
.PHONY: clean
