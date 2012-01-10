###########################################################
#ARCH	AIX  #serial dmpar
#

SFC=xlf_r 
SF90=xlf90  -qfree
SCC=xlc_r
FFLAGS=-O -g -qnosave -qarch=auto -qflttrap=zerodivide:enable -qsigtrap -C -qinitauto=FF911299 -q64 -C -qfullpath

DM_FC=mpxlf_r  
DM_F90=mpxlf_r  -qfree
DM_CC=mpxlc_r

FC              =       CONFIGURE_FC
F90             =       CONFIGURE_F90
CC              =       CONFIGURE_CC

CFLAGS=-DIBM4 -DVERBOSE -q64 -C -qfullpath

###########################################################
#ARCH	Linux i486 i586 i686, PGI compiler #serial dmpar
#
LDFLAGS         =      -Wl,-noinhibit-exec 

SFC=pgf90
SF90=pgf90 -Mfree -C 
SCC=pgcc
FFLAGS=-O0 -g -Kieee -pc 32 -Ktrap=fp -C -byteswapio $(LDFLAGS)
CPP             =      /lib/cpp -C -P

DM_FC=mpif90 -f90=pgf90
DM_F90=mpif90 -Mfree -f90=pgf90
DM_CC=mpicc 

FC              =       CONFIGURE_FC
F90             =       CONFIGURE_F90
CC              =       CONFIGURE_CC

CFLAGS=-O0 -DLINUX

###########################################################
#ARCH	Linux x86_64, PGI compiler # serial dmpar
#
LDFLAGS         =      -Wl,-noinhibit-exec 

SFC=pgf90
SF90=pgf90 -Mfree -C 
SCC=pgcc
FFLAGS=-O0 -g -Kieee -pc 64 -Ktrap=fp -C -byteswapio $(LDFLAGS)

DM_FC=mpif90 -f90=pgf90
DM_F90=mpif90 -Mfree -f90=pgf90
DM_CC=mpicc 

FC              =       CONFIGURE_FC
F90             =       CONFIGURE_F90
CC              =       CONFIGURE_CC

CFLAGS=-O0 -DLINUX


###########################################################
#ARCH	Linux i486 i586 i686, Intel compiler	# serial dmpar
#
LDFLAGS         =     -Wl,-noinhibit-exec 

SFC=ifort
SF90=ifort -free
SCC=icc
FFLAGS=-O3 -xT -fp-model precise -assume byterecl -convert big_endian -fpe0 -g -traceback $(LDFLAGS)

DM_FC=ifort
DM_F90=ifort -free
DM_CC=icc

FC              =       CONFIGURE_FC
F90             =       CONFIGURE_F90
CC              =       CONFIGURE_CC

CFLAGS=-O0 -DLINUX


###########################################################
#ARCH	Linux x86_64, Intel compiler	# serial dmpar
#
LDFLAGS         =     -Wl,-noinhibit-exec 

SFC=ifort
SF90=ifort -free
SCC=icc
FFLAGS=-O3 -xT -fp-model precise -assume byterecl -convert big_endian -fpe0 -g -traceback $(LDFLAGS)

DM_FC=ifort
DM_F90=ifort -free
DM_CC=icc

FC              =       CONFIGURE_FC
F90             =       CONFIGURE_F90
CC              =       CONFIGURE_CC

CFLAGS=-O0 -DLINUX

###########################################################
#ARCH	Linux i486 i586 i686, gfortran compiler # serial dmpar
#
LDFLAGS         =      -Wl,-noinhibit-exec

SFC=gfortran
SF90=gfortran -ffree-form
SCC=gcc
FFLAGS=-fconvert=big-endian -fno-second-underscore -fno-range-check -frecord-marker=4 $(LDFLAGS)


DM_FC=mpif90
DM_F90=mpif90 -ffree-form
DM_CC=mpicc

FC              =       CONFIGURE_FC
F90             =       CONFIGURE_F90
CC              =       CONFIGURE_CC

CFLAGS=-O0 -DLINUX -D_OPENMP

###########################################################
#ARCH	Linux x86_64, gfortran compiler # serial dmpar
#
LDFLAGS         =      -Wl,-noinhibit-exec

SFC=gfortran
SF90=gfortran -ffree-form
SCC=gcc
FFLAGS=-fconvert=big-endian -fno-second-underscore -fno-range-check $(LDFLAGS)

DM_FC=gfortran
DM_F90=gfortran -ffree-form
DM_CC=gcc

FC              =       CONFIGURE_FC
F90             =       CONFIGURE_F90
CC              =       CONFIGURE_CC

CFLAGS=-O0 -DLINUX -D_OPENMP


##################################################################
#ARCH	NULL
