#----------------------------------------------------------------------------
# makefile for the NCEP g2tmpl library installed within the UPP community 
# code distributed by the DTC
#----------------------------------------------------------------------------
SHELL = /bin/sh
LIB   = libg2tmpl.a
#
# configureation based on architecture and compiler set at UPP system
# directory level
include ../../../configure.upp

EXTRA_FFLAGS  = -c $(PROMOTION) -I$(INCMOD)
EXTRA_CFLAGS  = -c
EXTRA_ARFLAGS = 

OBJS_F    = grib2_all_tables_module.o gdtsec3.o \
            g2sec4_temp1.o g2sec4_temp11.o g2sec4_temp12.o
OBJS_C    = process_4dot2.o
OBJS      = $(OBJS_C) $(OBJS_F)

SRC_F     = $(OBJS_F:.o=.f)
SRC_C     = $(OBJS_C:.o=.c)

all: $(LIB)

$(LIB): $(OBJS)
	$(AR) $(ARFLAGS) $(EXTRA_ARFLAGS) $@ $(OBJS)
	$(CP) $(LIB) $(LIBDIR)
	$(CP) *.mod $(INCMOD)

clean:
	@for f in `ls *.mod` ; do \
		$(RM) $$f ; $(RM) $(INCMOD)/$$f ; \
	done
	$(RM) clib.h
	$(RM) $(LIBDIR)/$(LIB)
	$(RM) $(OBJS_F) $(OBJS_F77) $(OBJS_C) $(LIB)
#
# always build serial -- and this is new format code
.f.o:
	$(SF90) $(FFLAGS) $(EXTRA_FFLAGS) $<

.PHONY: clean
.IGNORE:
