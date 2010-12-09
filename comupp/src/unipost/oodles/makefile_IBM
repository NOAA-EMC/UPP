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
#         W3LIB    - w3lib
#         BACIO    - bacio lib
#         ESSL     - ESSL library
#         MASS     - MASS library
#         HPMLIB   - hpm lib
#         SEARCH   - library search location
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
#TARGET = ../exec/nceppost.x
TARGET = ncep_post
#
# CPP, Compiler, and Linker Options
#



WRFPATH = /nwprod/sorc/nam_nmm_real_fcst.fd


#NETCDFPATH = /usrx/local/netcdf.3.5.0
NETCDFPATH = /usrx/local/netcdf-4.0.1
FC       = mpxlf90_r
CPP      = /lib/cpp -P
ARCH     = auto
CPPFLAGS = 
#OPTS     = -O -qnosave -qarch=$(ARCH) -qmaxmem=-1 -NS2000
#OPTS     = -C -O -qnosave -qarch=$(ARCH) -qmaxmem=-1 -NS2000
#OPTS     = -O -g -qnosave -qarch=$(ARCH) -qmaxmem=-1 -NS2000 -bmaxdata:0x80000000
OPTS     = -O -g -qnosave -qarch=$(ARCH) -qmaxmem=-1 -NS2000
LIST     = 
FREE     = -qfree
#TRAPS    = -qflttrap=ov:und:zero:inv:inex -qcheck -qinitauto=FF
TRAPS    = 
PROFILE  = 
DEBUG    = -g
MEM      =
MAP      = -bmap:map -bloadmap:lm 
W3LIBDIR = /nwprod/lib
ESSL     = -lessl
MASS     = -lmass
#NCDLIBS = -L$(NETCDFPATH)/lib -lnetcdf
# new netcdf path
NCDLIBS = -L$(NETCDFPATH)/libsrc/.libs/ -lnetcdf
#NCDFFLAGS = -I$(NETCDFPATH)/include
# new netcdf path
NCDFFLAGS = -I$(NETCDFPATH)/libsrc

WRFFFLAGS = -I$(WRFPATH)/external/io_quilt
CRTMFFLAGS = -I/nwprod/lib/incmod/crtm2
W3FLAGS = -I/nwprod/lib/incmod/w3_4
SFCFLAG = -I/nwprod/lib/incmod/sfcio_4
WRFLIB    = $(WRFPATH)/main/libwrflib.a $(WRFPATH)/external/io_int/libwrfio_int.a $(WRFPATH)/external/io_netcdf/libwrfio_nf.a $(WRFPATH)/frame/pack_utils.o $(WRFPATH)/external/esmf_time_f90/libesmf_time.a $(WRFPATH)/external/RSL_LITE/librsl_lite.a
CRTMLIB = 

SEARCH   =
#
# Assemble Options
#
FFLAGS   = $(OPTS) $(LIST) $(TRAPS) $(PROFILE) $(DEBUG) $(NCDFFLAGS) $(WRFFLAGS) $(CRTMFFLAGS) $(W3FLAGS) $(SFCFLAG)
FFLAGST  = $(OPTS) $(LIST) $(FREE) $(TRAPS) $(PROFILE) $(DEBUG) $(NCDFFLAGS) $(WRFFLAGS) $(CRTMFFLAGS) $(W3FLAGS)
LDFLAGS  = $(MEM) $(MAP) $(SMP) $(PROFILE)
#LIBS     = $(ESSL) $(MASS) $(SEARCH) $(NCDLIBS) $(WRFLIB) -L$(W3LIBDIR) -lw3_4 -lbacio_4 -lsp_4 -lsigio_4 -lsfcio_4 -lcrtm2
LIBS     = $(ESSL) $(MASS) $(SEARCH) $(NCDLIBS) $(WRFLIB)\
           -L/nwpara/lib -lsigio_4 -lsfcio_4\
           -L/climate/save/wx20wa/gfsio/bacio -lbacio_4 -L/nwprod/lib -lsp_4 -lcrtm2 /global/save/wx20gg/bgrids/w3mods/w3/lib/libw3_4.a 
#
#
# Threaded object files
#
OBJST=	wrf_io_flags.o module_internal_header_util.o getVariable.o getIVariable.o getVariableB.o getIVariableN.o getVariableRSM.o \
	kinds_mod.o gfsio_module.o nemsio_module.o machine.o physcons.o \
	count_recs_wrf_binary_file.o inventory_wrf_binary_file.o \
	next_buf.o retrieve_index.o ZENSUN.o CLDFRAC_ZHAO.o \
	GFSPOST.o GETGBANDSCATTER.o 
#
# Non-threaded object files
#
OBJS=	VRBLS2D_mod.o VRBLS3D_mod.o MASKS_mod.o PMICRPH.o SOIL_mod.o \
        CMASSI.o CTLBLK.o GRIDSPEC.o LOOKUP.o PARAMR.o RHGRD.o RQSTFLD.o \
        cuparm.o params.o svptbl.o \
	BNDLYR.o  BOUND.o  CALCAPE.o  CALDWP.o  CALDRG.o CALHEL.o  CALLCL.o  \
	CALMCVG.o CALPOT.o  CALPW.o CALRH.o  CALRCH.o CALRH_GSD.o \
	CALSTRM.o CALTAU.o CALTHTE.o CALVIS.o CALVIS_GSD.o CALVOR.o CALWXT.o \
        CALWXT_RAMER.o CALWXT_BOURG.o CALWXT_REVISED.o \
        CALWXT_EXPLICIT.o CALWXT_DOMINANT.o \
	CLDRAD.o  CLMAX.o COLLECT.o  COLLECT_LOC.o \
	DEWPOINT.o \
	FDLVL.o FGAMMA.o FIXED.o  FRZLVL.o  FRZLVL2.o \
	GET_BITS.o  GRIBIT.o INITPOST.o LFMFLD.o  INITPOST_BIN.o \
	MAPSSLP.o MISCLN.o MIXLEN.o MDL2P.o MDLFLD.o MPI_FIRST.o  MPI_LAST.o \
	NGMFLD.o NGMSLP.o  OTLFT.o OTLIFT.o SLP_new.o SLP_NMM.o EXCH.o \
	PARA_RANGE.o              PROCESS.o INITPOST_NMM.o EXCH2.o \
	READCNTRL.o  SCLFLD.o  SERVER.o  SETUP_SERVERS.o SMOOTH.o SURFCE.o \
	SPLINE.o  TABLE.o  TABLEQ.o  TRPAUS.o  TTBLEX.o WETBULB.o WRFPOST.o \
        INITPOST_NMM_BIN.o CALMICT.o CALVIS.o MICROINIT.o GPVS.o MDL2SIGMA.o \
        ETCALC.o CANRES.o CALGUST.o WETFRZLVL.o SNFRAC.o MDL2AGL.o SNFRAC_GFS.o \
	INITPOST_RSM.o AVIATION.o DEALLOCATE.o INITPOST_NMM_BIN_MPIIO_IJK.o \
        CALPBL.o MDL2SIGMA2.o INITPOST_GFS.o CALRH_GFS.o LFMFLD_GFS.o CALRAD.o \
	CALRAD_WCLOUD.o MDL2THANDPV.o CALPBLREGIME.o POLEAVG.o \
	INITPOST_NEMS.o GETNEMSNDSCATTER.o ICAOHEIGHT.o INITPOST_GFS_NEMS.o 
#
# Includes
#
##INCLUDES= parm.tbl parmeta parmout parmsoil cuparm
#INCLUDES= parm.tbl cuparm
#
# Common Blocks
#
#COMMS=	LOOKUP.comm   RQSTFLD.comm   CTLBLK.comm  \
#        GRIDSPEC.comm CMASSI.comm RHGRD.comm

#DEPS= $(COMMS) $(INCLUDES)

.SUFFIXES:	.F .f .o

.F.f:
	$(CPP) $(CPPFLAGS) $< > $*.f

$(TARGET):	$(OBJST) $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJST) $(OBJS) $(LIBS)

$(OBJS):	$(DEPS)
	$(FC) $(FFLAGS) -c $<

$(OBJST):	$(DEPS)
	$(FC) $(FFLAGST) -c $<

clean:	
	/bin/rm -f  $(TARGET) *.lst *.o *.mod
#
