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

# Define the name of the executable


# build configuration determined before compile
include configure.wrf

COPY	= 	cp


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

WRF_COMPILE_START="`date`"

all:
	$(MAKE) wrf_guts
	@echo "wrf build started:   $(WRF_COMPILE_START)"
	@echo "wrf build completed:" `date`

wrf_guts: 
	@ echo '--------------------------------------'
	$(COPY) Registry/Registry.EM Registry/Registry
	( cd tools ; $(MAKE) standard.exe )
	(            $(MAKE) io_only )
	(            $(MAKE) MODULE_DIRS="$(ALL_MODULES)" toolsdir )
	( cd frame ; $(MAKE) module_driver_constants.o pack_utils.o module_machine.o module_internal_header_util.o wrf_debug.o )
	( cd frame ; $(AR) $(ARFLAGS) ../main/libwrflib.a module_driver_constants.o pack_utils.o module_machine.o  \
		module_internal_header_util.o module_wrf_error.o wrf_debug.o )

io_only:  esmf_time wrfio_nf   \
          wrf_ioapi_includes wrfio_grib_share wrfio_grib1 wrfio_int fftpack

esmf_time :
	( cd external/esmf_time_f90 ; \
          make $(J) FC="$(SFC) $(PROMOTION) $(FCDEBUG) $(FCBASEOPTS)" RANLIB="$(RANLIB)" \
          CPP="$(CPP) -Iinc -I. $(ARCHFLAGS) $(TRADFLAG)" AR="$(AR)" ARFLAGS="$(ARFLAGS)" )

wrfio_nf :
	( cd external/io_netcdf ; \
          make $(J) NETCDFPATH="$(NETCDFPATH)" RANLIB="$(RANLIB)" CPP="$(CPP)" \
          CC="$(SCC)" CFLAGS="$(CFLAGS)" \
          FC="$(SFC) $(PROMOTION) $(OMP) $(FCFLAGS)" TRADFLAG="$(TRADFLAG)" AR="$(AR)" ARFLAGS="$(ARFLAGS)" )


wrf_ioapi_includes :
	( cd external/ioapi_share ; \
          $(MAKE) NATIVE_RWORDSIZE="$(NATIVE_RWORDSIZE)" RWORDSIZE="$(RWORDSIZE)" AR="$(AR)" ARFLAGS="$(ARFLAGS)" )

wrfio_grib_share :
	( cd external/io_grib_share ; \
          make $(J) CC="$(SCC)" CFLAGS="$(CFLAGS)" RM="$(RM)" RANLIB="$(RANLIB)" CPP="$(CPP)" \
          FC="$(SFC) $(PROMOTION) -I. $(FCDEBUG) $(FCBASEOPTS) $(FCSUFFIX)" TRADFLAG="$(TRADFLAG)" AR="$(AR)" ARFLAGS="$(ARFLAGS)" archive)

wrfio_grib1 :
	( cd external/io_grib1 ; \
          make $(J) CC="$(SCC)" CFLAGS="$(CFLAGS)" RM="$(RM)" RANLIB="$(RANLIB)" CPP="$(CPP)" \
          FC="$(SFC) $(PROMOTION) -I. $(FCDEBUG) $(FCBASEOPTS) $(FCSUFFIX)" TRADFLAG="$(TRADFLAG)" AR="$(AR)" ARFLAGS="$(ARFLAGS)" archive)

wrfio_int :
	( cd external/io_int ; \
          make $(J) CC="$(CC)" CFLAGS_LOCAL="$(CFLAGS_LOCAL)" RM="$(RM)" RANLIB="$(RANLIB)" CPP="$(CPP)" \
          FC="$(FC) $(PROMOTION) $(FCDEBUG) $(FCBASEOPTS) $(OMP)" FGREP="$(FGREP)" \
          TRADFLAG="$(TRADFLAG)" AR="$(AR)" ARFLAGS="$(ARFLAGS)" ARCHFLAGS="$(ARCHFLAGS)" all )
fftpack :
	( cd external/fftpack/fftpack5 ; \
          make $(J) FC="$(SFC)" FFLAGS="$(PROMOTION) $(FCDEBUG) $(FCBASEOPTS)" RANLIB="$(RANLIB)" AR="$(AR)" \
          ARFLAGS="$(ARFLAGS)" CPP="$(CPP)" CPPFLAGS="$(CPPFLAGS)" RM="$(RM)" )

toolsdir :
	  ( cd tools ; $(MAKE) CC_TOOLS_CFLAGS="$(CC_TOOLS_CFLAGS)" CC_TOOLS="$(CC_TOOLS) -DIWORDSIZE=$(IWORDSIZE) -DMAX_HISTORY=$(MAX_HISTORY)" ) ; \

