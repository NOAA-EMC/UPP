SHELL = /bin/sh

################################################################################
# 
#     Makefile for NCEP Post
#
#     Use:
#     make         -  build the executable
#     make clean   -  start with a clean slate
#
#     The following macros will be of interest:
#
#         TARGET   - name of the executable
#         FC       - name of Fortran compiler
#         CPP      - name of CPP
#         ARCH     - architecture
#         CPPFLAGS - CPP flags
#         OPTS     - compiler code optimizations
#         LIST     - source listing
#         SMP      - threading
#         TRAPS    - runtime traps for floating point exceptions
#         PROFILE  - source code profiling ( -pg )
#         DEBUG    - -g
#         MEM      - user data area and stack size
#         MAP      - load map
# 
#         This makefile was created based off the NCEP unipost makefile
#           changes were made to have in conform to the architecture
#           supprted by the DTC community code
#
#         This version for eta_post with more intelligent memory allocation
#         Jim Tuccillo   Feb 2001
# 
#         This version for eta_post with asynchronous I/O server.   
#         Jim Tuccillo   June 2001
#
#################################################################################
#
# Define the name of the executable
#
TARGET = ncep_post

#
# build configuration determined before compile
include ../../configure.upp

#
# Install directory for executable
BINDIR     = ../../bin
#
# directories for shared resources
INCMOD      = ../../include
LOCALINC    = -I$(INCMOD) -I$(INCMOD)/upp_crtm
NCDFINC     = -I$(NETCDFPATH)/include
WRFINC      = -I$(WRF_DIR)/external/io_quilt -I$(WRF_DIR)/frame

LIBDIR      = -L../../lib
UPPLIBS     = -lbacio -lupp_crtm -lsigio -lsfcio -lsp -lmersenne -lw3
NCDFLIBS    = -L$(NETCDFPATH)/lib $(NETCDFLIBS)
WRFLIB      = -L$(WRF_DIR)/main -lwrflib                          \
              -L$(WRF_DIR)/external/io_int -lwrfio_int            \
              -L$(WRF_DIR)/external/io_netcdf -lwrfio_nf          \
              -L$(WRF_DIR)/external/esmf_time_f90 -lesmf_time     \
              -L$(WRF_DIR)/external/RSL_LITE  -lrsl_lite          \
              -L$(WRF_DIR)/external/io_grib1 -lio_grib1           \
              -L$(WRF_DIR)/external/io_grib_share -lio_grib_share \
              -L$(WRF_DIR)/external/fftpack/fftpack5 -lfftpack
LIBS        = $(LIBDIR) $(UPPLIBS) $(NCDFLIBS) $(WRFLIB)

MODULES     = ../NCEP_modules/kinds_mod.o ../NCEP_modules/constants_mod.o \
               $(WRF_DIR)/frame/module_internal_header_util.o             \
               $(WRF_DIR)/frame/pack_utils.o

#
# Compilation / Link Flag Configuration
EXTRA_CPPFLAGS = -DLINUX
EXTRA_FFLAGS   = -c $(LOCALINC) $(NETCDFINC) $(WRFINC)
EXTRA_LDFLAGS  = $(LIBS)

# -----------
# Threaded object files
# -----------
OBJS_FT = wrf_io_flags.o getVariable.o getIVariable.o getVariableB.o getIVariableN.o getVariableRSM.o \
          gfsio_module.o nemsio_module.o machine.o physcons.o count_recs_wrf_binary_file.o \
          inventory_wrf_binary_file.o next_buf.o retrieve_index.o ZENSUN.o CLDFRAC_ZHAO.o GFSPOST.o \
          GETGBANDSCATTER.o rsearch.o

# -----------
# Non-threaded object files
# -----------
OBJS_F =	 VRBLS2D_mod.o VRBLS3D_mod.o MASKS_mod.o PMICRPH.o SOIL_mod.o CMASSI.o CTLBLK.o GRIDSPEC.o \
          LOOKUP.o PARAMR.o RHGRD.o RQSTFLD.o cuparm.o params.o svptbl.o BNDLYR.o  BOUND.o  CALCAPE.o \
          CALDWP.o  CALDRG.o CALHEL.o  CALLCL.o  CALMCVG.o CALPOT.o  CALPW.o CALRH.o  CALRCH.o \
          CALRH_GSD.o CALSTRM.o CALTAU.o CALTHTE.o CALVIS.o CALVIS_GSD.o CALVOR.o CALWXT.o \
          CALWXT_RAMER.o CALWXT_BOURG.o CALWXT_REVISED.o CALWXT_EXPLICIT.o CALWXT_DOMINANT.o CLDRAD.o \
          CLMAX.o COLLECT.o  COLLECT_LOC.o DEWPOINT.o FDLVL.o  FGAMMA.o FIXED.o  FRZLVL.o  FRZLVL2.o \
          GET_BITS.o  GRIBIT.o INITPOST.o LFMFLD.o  INITPOST_BIN.o MAPSSLP.o MISCLN.o MIXLEN.o MDL2P.o \
          MDLFLD.o MPI_FIRST.o  MPI_LAST.o NGMFLD.o NGMSLP.o  OTLFT.o OTLIFT.o SLP_new.o SLP_NMM.o \
          EXCH.o PARA_RANGE.o PROCESS.o INITPOST_NMM.o EXCH2.o READCNTRL.o  SCLFLD.o  SERVER.o \
          SETUP_SERVERS.o SMOOTH.o SURFCE.o SPLINE.o  TABLE.o  TABLEQ.o  TRPAUS.o  TTBLEX.o WETBULB.o \
          WRFPOST.o INITPOST_NMM_BIN.o CALMICT.o MICROINIT.o GPVS.o MDL2SIGMA.o ETCALC.o CANRES.o \
          CALGUST.o WETFRZLVL.o SNFRAC.o MDL2AGL.o SNFRAC_GFS.o INITPOST_RSM.o AVIATION.o DEALLOCATE.o \
          INITPOST_NMM_BIN_MPIIO_IJK.o CALPBL.o MDL2SIGMA2.o INITPOST_GFS.o CALRH_GFS.o LFMFLD_GFS.o \
          CALRAD.o CALRAD_WCLOUD.o MDL2THANDPV.o CALPBLREGIME.o POLEAVG.o INITPOST_NEMS.o \
          GETNEMSNDSCATTER.o ICAOHEIGHT.o INITPOST_GFS_NEMS.o 

OBJS   = $(OBJS_FT) $(OBJS_F)

# -----------
# Targets
# -----------
all: $(TARGET)

$(TARGET):	$(OBJS_F)
	$(FC) -o $@ $(LDFLAGS) $(EXTRA_LDFLAGS) $(OBJS) $(MODULES)
	$(CP) $@ $(BINDIR)

#
# This insures a dependency found in some files -- watch file order above remains -- should
# be done w/ dependencies
$(OBJS_F): $(OBJS_FT)

clean:	
	@echo -e "\n<><><><> CLEAN <><><><>\n$@ in `pwd`"
	$(RM) $(TARGET) $(OBJS) *.lst *.mod
	@for f in `ls -1 *.F|sed "s/.F$$/.f/"` ; do \
       $(RM) $$f   ; \
   done

distclean: clean

.SUFFIXES:
.SUFFIXES:	.F .f .o

.F.o:
	$(CPP) $(CPP_FLAGS) $< > $*.f
	$(F90)  -c $(FFLAGS) $(EXTRA_FFLAGS) $*.f

.f.o: 
	$(F90)  -c $(FFLAGS) $(EXTRA_FFLAGS) $<
