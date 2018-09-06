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

#Make objects for io_int for binary i/o
INTOBJS = io_int.o io_int_idx.o module_io_int_idx.o pack_utils.o

M4      = m4 -Uinclude -Uindex -Ulen


# directories for shared resources
LLIBDIR     = -L$(LIBDIR)
WRFEXTLIBS  = $(WRF_LIB) $(WRF_LIB_EXTRA) $(WRF_LIB2)
NCDFLIBS    = -L$(NETCDFPATH)/lib $(NETCDFLIBS)

LIBS        = $(LLIBDIR) $(WRFEXTLIBS) $(NCDFLIBS)

MODULES     = $(WRF_MODS)

#In WRF configure.defaults this is set separately, but for all of UPP's supported platforms they are the same
CC_TOOLS     = $(SCC)

#UPP has a lot of different flag/name conventions from WRF, somehow
FCFLAGS      = $(FOPT) $(FFLAGS)
FCDEBUG      = # $(FDEBUG) # Uncomment this for debug compile

#Hard-coding this for now, but may be platform-dependent
RWORDSIZE = $(NATIVE_RWORDSIZE)

#Aliases specific to WRF_IO_INT
FCINT=$(F90) $(FCDEBUG) $(OMP)
CPPINT    = $(CPP) -P -nostdinc -traditional-cpp $(TRADFLAG) $(ARCHFLAGS)


#Other settings for this makefile
WRF_COMPILE_START="`date`"

#Makefile rules
.SUFFIXES: .F90 .F .f .f90 .c .h .o .code

all:
	$(MAKE) wrf_guts
	@echo "wrf build started:   $(WRF_COMPILE_START)"
	@echo "wrf build completed:" `date`

wrf_guts: 
	@ echo '--------------------------------------'
	(            $(MAKE) io_only )

io_only:  wrfio_nf   \
          libwrfio_int.a

wrfio_nf :
	( cd external/io_netcdf ; \
          make $(J) NETCDFPATH="$(NETCDFPATH)" RANLIB="$(RANLIB)" CPP="$(CPP) -P -nostdinc -traditional-cpp" \
          CC="$(SCC)" CFLAGS="$(CFLAGS)" \
          FC="$(SF90) $(OMP) $(FCFLAGS)" TRADFLAG="$(TRADFLAG)" AR="$(AR)" ARFLAGS="$(ARFLAGS)" )


#### IO_INT section ####

libwrfio_int.a: $(INTOBJS)
	$(RM) libwrfio_int.a
	$(AR) $(ARFLAGS) libwrfio_int.a $(INTOBJS)
	$(RANLIB) libwrfio_int.a

io_int.f: io_int.F90 module_internal_header_util.o
	$(CPPINT) $*.F90 | $(M4) - > $@

io_int.o: io_int.f intio_tags.h
	$(FCINT) $(FCFLAGS) -o $@ -c $*.f

io_int_idx_tags.h: intio_tags.h
	awk '{print "#define", toupper($$4), $$6}' < intio_tags.h > $@

io_int_idx.o: io_int_idx.c io_int_idx.h io_int_idx_tags.h
	$(CC) -o $@ -c $(CFLAGS_LOCAL) $*.c

module_io_int_idx.o: module_io_int_idx.f
	$(FCINT) $(FCFLAGS) -o $@ -c $*.f

module_io_int_idx.f: module_io_int_idx.F90
	$(CPPINT) $*.F90 > $@

module_internal_header_util.f: module_internal_header_util.F pack_utils.o
	-$(RM) $@
	$(CPPINT) module_internal_header_util.F > $@

module_internal_header_util.o: module_internal_header_util.f
	$(FCINT) $(FCFLAGS) -o $@ -c $*.f

pack_utils.o :
	$(CC) -c $(CFLAGS) pack_utils.c


#### Other makerules ####

clean : 
	  ( cd external ; make -i superclean )

