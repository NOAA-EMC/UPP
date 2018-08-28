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

LLIBDIR     = -L$(LIBDIR)
WRFEXTLIBS  = $(WRF_LIB) $(WRF_LIB_EXTRA) $(WRF_LIB2)
NCDFLIBS    = -L$(NETCDFPATH)/lib $(NETCDFLIBS)

LIBS        = $(LLIBDIR) $(WRFEXTLIBS) $(NCDFLIBS)

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
	(            $(MAKE) io_only )
	( cd frame ; $(MAKE) pack_utils.o \
                     FC="$(SF90) $(FCFLAGS)" CPP="$(CPP) $(CPP_WRFIO)" SFC="$(SF90) $(FCFLAGS)" CPP="$(CPP) $(CPP_WRFIO)")
	( cd frame ; $(AR) $(ARFLAGS) ../libwrflib.a pack_utils.o )

io_only:  wrfio_nf   \
          wrfio_int

wrfio_nf :
	( cd external/io_netcdf ; \
          make $(J) NETCDFPATH="$(NETCDFPATH)" RANLIB="$(RANLIB)" CPP="$(CPP) $(CPP_WRFIO)" \
          CC="$(SCC)" CFLAGS="$(CFLAGS)" \
          FC="$(SF90) $(OMP) $(FCFLAGS)" TRADFLAG="$(TRADFLAG)" AR="$(AR)" ARFLAGS="$(ARFLAGS)" )


wrfio_int :
	( cd external/io_int ; \
          make $(J) CC="$(CC)" CFLAGS_LOCAL="$(CFLAGS_LOCAL)" RM="$(RM)" RANLIB="$(RANLIB)" CPP="$(CPP) $(CPP_WRFIO)" \
          FC="$(F90) $(FCDEBUG) $(FCBASEOPTS) $(OMP)" FGREP="$(FGREP)" \
          TRADFLAG="$(TRADFLAG)" AR="$(AR)" ARFLAGS="$(ARFLAGS)" ARCHFLAGS="$(ARCHFLAGS)" all )

clean : 
	  ( cd external ; make -i superclean )

