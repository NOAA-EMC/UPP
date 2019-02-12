################################################################################
#
#     Makefile for external WRF code. This should be built first, as it is 
#     needed by other parts of the UPP build
#
#     Use:
#     make         -  build the executable
#     make clean   -  start with a clean slate
#
#################################################################################

# Include configuration info
include ../../../configure.upp
IWORDSIZE = 4
LWORDSIZE = 4

##############################################################################

#### Native size (in bytes) of Fortran REAL data type on this architecture ####
#### Note:  to change real wordsize (for example, to promote REALs from.
####        4-byte to 8-byte), modify the definition of RWORDSIZE in the.
####        section following "Architecture specific settings".  Do not.
####        change NATIVE_RWORDSIZE as is it architecture-specific...
NATIVE_RWORDSIZE = 4
RWORDSIZE       =       $(NATIVE_RWORDSIZE)

#### NETCDF4 pieces

NETCDF4_IO_OPTS = -DUSE_NETCDF4_FEATURES -DWRFIO_NCD_LARGE_FILE_SUPPORT

#Make objects for netCDF i/o
NCOBJS = wrf_io.o field_routines.o

M4      = m4 -Uinclude -Uindex -Ulen


# directories for shared resources
WRFEXTLIBS  = $(WRF_LIB) $(WRF_LIB_EXTRA) $(WRF_LIB2)
NCDFLIBS    = -L$(NETCDFPATH)/lib $(NETCDFLIBS)

MODULES     = $(WRF_MODS)

#UPP has a lot of different flag/name conventions from WRF, somehow
FCFLAGS      = $(FOPT) $(FFLAGS)
FCDEBUG      = # $(FDEBUG) # Uncomment this for debug compile

#Hard-coding this for now, but may be platform-dependent
RWORDSIZE = $(NATIVE_RWORDSIZE)

#Aliases specific to WRF_IO_INT
FCINT  = $(F90) $(FCDEBUG) $(OMP)
CPPINT = $(CPP) -P -nostdinc -traditional-cpp $(TRADFLAG) $(ARCHFLAGS)

#Aliases specific to WRF_IO_NC
FCNC   = $(SF90) $(OMP) $(FCFLAGS)
CPPNC  = $(CPP) -P -nostdinc -traditional-cpp $(TRADFLAG)
CODE    = ext_ncd_get_dom_ti.code ext_ncd_get_var_td.code ext_ncd_get_var_ti.code ext_ncd_put_dom_ti.code ext_ncd_put_var_td.code ext_ncd_put_var_ti.code transpose.code

#Other settings for this makefile
WRF_COMPILE_START="`date`"

#Makefile rules
.SUFFIXES: .F90 .F .f .f90 .c .h .o .code

all:
	$(MAKE) io_only
	$(CP) *.a $(LIBDIR)
	$(CP) *.mod $(INCMOD)
	@echo "wrf build started:   $(WRF_COMPILE_START)"
	@echo "wrf build completed:" `date`

io_only:  libwrfio_nf.a

#### IO_NETCDF section ####

libwrfio_nf.a:   $(NCOBJS) $(CODE)
		/bin/rm -f $@
		$(AR) cr libwrfio_nf.a $(NCOBJS) ; \
		$(RANLIB) $@

wrf_io.o:       wrf_io.F90 $(CODE)
		grep nf_format_64bit $(NETCDFPATH)/include/netcdf.inc ;\
                a=$$? ; export a ; \
		if [ $$a -a "$$WRFIO_NCD_LARGE_FILE_SUPPORT" = "1" ] ; then \
                  $(CPPNC) -DWRFIO_NCD_LARGE_FILE_SUPPORT wrf_io.F90 | $(M4) - > wrf_io.f ; \
                else \
                  $(CPPNC) wrf_io.F90 | $(M4) - > wrf_io.f ; \
                fi
		$(FCNC) -w -o $@ $(FCFLAGS) -I$(NETCDFPATH)/include -c wrf_io.f

field_routines.o: field_routines.F90 wrf_io.o
		$(CPPNC) field_routines.F90 > field_routines.f
		$(FCNC) -w -o $@ $(FCFLAGS) -I$(NETCDFPATH)/include -c field_routines.f


#### Other makerules ####

clean : 
	( $(RM) *.f *.o *.obj *.i *.mod libwrfio_nf.a )

.f90.o:
	$(RM) $@
	$(FC) -o $@ -c $(FCFLAGS) $(PROMOTION) $(FCSUFFIX) $*.f90

.c.o:
	$(RM) $@
	$(CC) -o $@ -c $(CFLAGS) $*.c

