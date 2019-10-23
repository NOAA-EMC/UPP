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
# Object files for this library - all are f90
OBJS_F   =
OBJS_F77 =                                    \
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

clean:
	$(RM) $(LIBDIR)/$(LIB)
	$(RM) $(OBJS) $(LIB)

# Always serial build
.f.o:
	$(SFC) $(FFLAGS) $(EXTRA_FFLAGS) $<

.IGNORE:
.PHONY: clean
