SHELL = /bin/sh

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
include ../configure.upp
include wrf_makerules

# directories for shared resources
LOCALINC    = -I$(INCMOD) -I$(INCMOD)/crtm2
NCDFINC     = -I$(NETCDFPATH)/include
WRFINC      = -I$(WRF_MODS)
GRIB2INC    = -I$(GRIB2SUPT_INC)

LLIBDIR     = -L$(LIBDIR)
GRIB2LIBS   = -lg2 -lg2tmpl -lxmlparse $(GRIB2SUPT_LIB)
UPPLIBS     = -lnemsio -lsigio -lsfcio -lgfsio -lsp -lw3nco -lw3emc -lbacio -lCRTM $(SERIAL_MPI_LIB)
WRFEXTLIBS  = $(WRF_LIB) $(WRF_LIB_EXTRA) $(WRF_LIB2)
NCDFLIBS    = -L$(NETCDFPATH)/lib $(NETCDFLIBS)

LIBS        = $(LLIBDIR) $(UPPLIBS) $(GRIB2LIBS) $(WRFEXTLIBS) $(NCDFLIBS)

MODULES     = $(WRF_MODS)

ALL_MODULES =                           \
               $(EM_MODULE_DIR)         \
               $(NMM_MODULES)           \
               $(EXP_MODULES)           \
               $(INCLUDE_MODULES)

#These CPP flags (not used by UPP code) are needed for some WRFIO stuff
CPP_WRFIO = -P -nostdinc -traditional-cpp

#In WRF configure.defaults this is set separately, but for all of UPP's supported platforms they are the same
CC_TOOLS     = $(SCC)

#UPP has a lot of different flag/name conventions from WRF, somehow
FCFLAGS      = $(FOPT) $(FFLAGS)
FCDEBUG      = # $(FDEBUG) # Uncomment this for debug compile

#Not sure why this is necessary for the toolsdir build, but it is
CC_TOOLS_CFLAGS = -DNMM_CORE=$(WRF_NMM_CORE)

#Hard-coding this for now, but may be platform-dependent
RWORDSIZE = $(NATIVE_RWORDSIZE)

# Should include this macro in configure.upp
FGREP = fgrep -iq

#Other settings for this makefile
WRF_COMPILE_START="`date`"


#Makefile rules

all:
	$(MAKE) wrf_guts
	@echo "wrf build started:   $(WRF_COMPILE_START)"
	@echo "wrf build completed:" `date`

wrf_guts: 
	@ echo '--------------------------------------'
	( cd tools ; $(MAKE) standard.exe )
	(            $(MAKE) io_only )
	( cd frame ; $(MAKE) module_driver_constants.o pack_utils.o module_machine.o module_internal_header_util.o wrf_debug.o \
                     FC="$(SF90) $(FCFLAGS)" CPP="$(CPP) $(CPP_WRFIO)" SFC="$(SF90) $(FCFLAGS)" CPP="$(CPP) $(CPP_WRFIO)")
	( cd frame ; $(AR) $(ARFLAGS) ../libwrflib.a module_driver_constants.o pack_utils.o module_machine.o  \
		module_internal_header_util.o module_wrf_error.o wrf_debug.o )

io_only:  esmf_time wrfio_nf   \
          wrf_ioapi_includes wrfio_grib_share wrfio_grib1 wrfio_int fftpack

esmf_time :
	( cd external/esmf_time_f90 ; \
          make $(J) FC="$(SF90) $(FCBASEOPTS)" RANLIB="$(RANLIB)" \
          CPP="$(CPP) $(CPP_WRFIO) -I. $(ARCHFLAGS) $(TRADFLAG)" AR="$(AR)" ARFLAGS="$(ARFLAGS)" )

wrfio_nf :
	( cd external/io_netcdf ; \
          make $(J) NETCDFPATH="$(NETCDFPATH)" RANLIB="$(RANLIB)" CPP="$(CPP) $(CPP_WRFIO)" \
          CC="$(SCC)" CFLAGS="$(CFLAGS)" \
          FC="$(SF90) $(OMP) $(FCFLAGS)" TRADFLAG="$(TRADFLAG)" AR="$(AR)" ARFLAGS="$(ARFLAGS)" )


wrf_ioapi_includes :
	( cd external/ioapi_share ; \
          $(MAKE) NATIVE_RWORDSIZE="$(NATIVE_RWORDSIZE)" RWORDSIZE="$(RWORDSIZE)" AR="$(AR)" ARFLAGS="$(ARFLAGS)" )

wrfio_grib_share :
	( cd external/io_grib_share ; \
          make $(J) CC="$(SCC)" CFLAGS="$(CFLAGS)" RM="$(RM)" RANLIB="$(RANLIB)" CPP="$(CPP) $(CPP_WRFIO)" \
          FC="$(SFC) -I. $(FCDEBUG) $(FCBASEOPTS) $(FCSUFFIX)" TRADFLAG="$(TRADFLAG)" AR="$(AR)" ARFLAGS="$(ARFLAGS)" archive)

wrfio_grib1 :
	( cd external/io_grib1 ; \
          make $(J) CC="$(SCC)" CFLAGS="$(CFLAGS)" RM="$(RM)" RANLIB="$(RANLIB)" CPP="$(CPP) $(CPP_WRFIO)" \
          FC="$(SFC) -I. $(FCDEBUG) $(FCBASEOPTS) $(FCSUFFIX)" TRADFLAG="$(TRADFLAG)" AR="$(AR)" ARFLAGS="$(ARFLAGS)" archive)

wrfio_int :
	( cd external/io_int ; \
          make $(J) CC="$(CC)" CFLAGS_LOCAL="$(CFLAGS_LOCAL)" RM="$(RM)" RANLIB="$(RANLIB)" CPP="$(CPP) $(CPP_WRFIO)" \
          FC="$(F90) $(FCDEBUG) $(FCBASEOPTS) $(OMP)" FGREP="$(FGREP)" \
          TRADFLAG="$(TRADFLAG)" AR="$(AR)" ARFLAGS="$(ARFLAGS)" ARCHFLAGS="$(ARCHFLAGS)" all )
fftpack :
	( cd external/fftpack/fftpack5 ; \
          make $(J) FC="$(SFC)" FFLAGS="$(FCDEBUG) $(FCBASEOPTS)" RANLIB="$(RANLIB)" AR="$(AR)" \
          ARFLAGS="$(ARFLAGS)" CPP="$(CPP) $(CPP_WRFIO)" CPPFLAGS="$(CPPFLAGS)" RM="$(RM)" )

clean : 
	  ( cd external ; make -i superclean )
	  ( cd external/io_grib1/WGRIB ; make clean )

