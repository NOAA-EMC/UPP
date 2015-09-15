#!/bin/ksh
set -x
mac=$(hostname | cut -c1-1)
mac2=$(hostname | cut -c1-2)
################################# options ###############################################
3export CLEAN=NO                                 # uncomment this if you don't want to clean
                                                 # before compiling
#debug=YES                                       # to turn on debug mode - defaults to NO
 make_post_lib=YES                               # to create post library - defaults to NO
 make_post_exec=YES                              # to create post executable - defaults to YES
 make_nowrf=YES                                  # to compile with wrf stub instead of WRF lib
#BMPYXML=_bmpyxml                                # to use original bumpy xml file
                                                 # make sure to clean when changing thisi                                                 # variable BMPXML 
################################# options ###############################################
#
if [ $mac2 = ga ] ; then                         # For GAEA
 machine=gaea
 center=${center:-ncep}
elif [ $mac2 = tf ] ; then                       # For Theia
 machine=theia
elif [ $mac = z -o $mac = h -o $mac = f ] ; then # For ZEUS
 machine=zeus
elif [ $mac = t -o $mac = e -o $mac = g ] ; then # For WCOSS
 machine=wcoss
fi
debug=${debug:-NO}
make_post_lib=${make_post_lib:-NO}
make_post_exec=${make_post_exec:-YES}
BMPYXML=${BMPYXML:-""}
if [ $machine = wcoss ] ; then
  export NETCDFPATH="/usrx/local/NetCDF/3.6.3"
  export WRFPATH="/nwprod/sorc/wrf_shared.v1.1.0"
  export NWPROD="/nwprod"
  export XMLPATH=$NWPROD
  export IPPATH=$NWPROD
  export SPPATH=/usrx/local/nceplibs
  ecport BACIOPATH=/usrx/local/nceplibs
  export ipv=""
  export spv=_v2.0.2p
  export crtmv=2.0.6
  export crtmv_inc=$crtmv
  export xmlv=_v2.0.0
  export baciov=_v2.0.1p
  export FC=mpiifort
  export CPP="/lib/cpp -P"
  export CPPFLAGS="-DLINUX"
  export CC=cc
  if [ $debug = YES ] ; then
    export OPTS="-O0 -openmp "
#   export DEBUG="-g -check all -ftrapuv -convert big_endian -fp-stack-check -fstack-protector -heap-arrays -recursive -traceback"
    export DEBUG="-g -traceback -convert big_endian -ftrapuv -check bounds -check format -check output_conversion -check pointers -check uninit -fp-stack-check"
  else
    export OPTS="-O3 -convert big_endian -fp-model source -openmp -xAVX"
    export DEBUG=""
  fi
  export LIST=""
  export FREE="-FR"
  export TRAPS=""
  export PROFILE=""
elif [ $machine = zeus ] ; then
  export NETCDFPATH="/apps/netcdf/3.6.3/intel"
  export WRFPATH="/scratch2/portfolios/NCEPDEV/meso/save/Dusan.Jovic/WRFV3"
  export NWPROD="/contrib/nceplibs/nwprod"
  export XMLPATH="/home/Hui-Ya.Chuang"
  export IPPATH=$NWPROD
  export SPPATH=$NWPROD
  export ipv=""
  export spv=_v2.0.1
  export crtmv=2.0.7
  export FC="ifort -lmpi"
  export CPP="/lib/cpp -P"
  export CC=cc
  export ARCH=""
  export CPPFLAGS="-DLINUX"
  if [ $debug = YES ] ; then
    export OPTS="-O0 -openmp -g"
    export DEBUG="-g -check all -ftrapuv -convert big_endian -fp-stack-check -fstack-protector -heap-arrays -recursive -traceback"
  else
    export export OPTS="-O3 -convert big_endian -traceback -g -fp-model source -openmp"
    export DEBUG=""
  fi
  export LIST=""
  export FREE="-FR"
  export TRAPS=""
  export PROFILE=""
elif [ $machine = theia ] ; then
  export NETCDFPATH="/apps/netcdf/4.3.0-intel"
# export WRFPATH="/scratch4/NCEPDEV/meso/save/Dusan.Jovic/WRFV3"
  export WRFPATH="/scratch4/NCEPDEV/global/save/Shrinivas.Moorthi/theia/nceplibs/nwprod/lib/sorc/WRFV3"

  export NWPROD="/scratch4/NCEPDEV/global/save/Shrinivas.Moorthi/theia/nceplibs/nwprod"
  export ipv=_v2.0.3
  export spv=""
##export spv=_v2.0.1
  export crtmv=2.0.7
  export gfsiov=""
  export w3ev=_v2.1.0
  export w3nv=""
  export xmlv=_v2.0.0
  export g2tv=""
  export baciov=_v2.1.0

# export NWPROD="/scratch3/NCEPDEV/nwprod"
# export ipv=_v2.0.0
# export spv=_v2.0.2
# export crtmv=2.0.6
# export crtmv=2.1.3
# export gfsiov=_v1.1.0
# export w3ev=_v2.0.5
# export w3nv=_v2.0.6
# export xmlv=_v2.0.0
# export g2tv=_v1.3.0
# export baciov=_v2.0.1

  export XMLPATH=$NWPROD
  export IPPATH=$NWPROD
  export SPPATH=$NWPROD
  export BACIOPATH=$NWPROD/lib

  export FC=mpiifort
  export CPP="/lib/cpp -P"
  export CC=cc
  export ARCH=""
  export CPPFLAGS="-DLINUX"
  if [ $debug = YES ] ; then
    export OPTS="-O0 -openmp -g"
    export DEBUG="-g -check all -ftrapuv -convert big_endian -fp-stack-check -fstack-protector -heap-arrays -recursive -traceback"
  else
    export export OPTS="-O3 -convert big_endian -traceback -g -fp-model source -openmp"
#   export export OPTS="-O2 -convert big_endian -traceback -g -fp-model source "
#   export export OPTS="-O2 -convert big_endian -traceback -g -fp-model source -openmp"
    export DEBUG=""
  fi
  export LIST=""
  export FREE="-FR"
  export TRAPS=""
  export PROFILE=""
  export make_nowrf=${make_nowrf:-YES}
elif [ $machine = gaea ] ; then
# export NETCDFPATH="/opt/cray/netcdf/4.1.1.0/netcdf-intel"
  export NETCDFPATH="/opt/cray/netcdf/4.2.0/intel/120/"
  export WRFPATH="/lustre/f1/unswept/ncep/Shrinivas.Moorthi/nceplibs/nwprod/lib/sorc/WRFV3"
# export WRFPATH="/lustre/f1/unswept/ncep/Shrinivas.Moorthi/nceplibs/nwprod/lib/sorc/wrf_shared.v1.1.0"
  export NWPROD="/lustre/f1/unswept/ncep/Shrinivas.Moorthi/nceplibs/nwprod"
  export IPPATH=$NWPROD
  export SPPATH=$NWPROD
  export ipv=""
  export spv=_v2.0.1
  export xmlv=_v2.0.0
  export FC=ftn
  export CPP="/lib/cpp -P"
  export ARCH=""
  export CPPFLAGS="-DLINUX"
  export CC=icc
  if [ $debug = YES ] ; then
    export OPTS="-O0 -g"
    export DEBUG="-g -check all -ftrapuv -convert big_endian -fp-stack-check -fstack-protector -heap-arrays -recursive -traceback"
  else
    export export OPTS="-O3 -convert big_endian -traceback -g -fp-model source"
    export DEBUG=""
  fi
  export LIST=""
  export FREE=-FR
  export TRAPS=""
  export PROFILE=""

  export gfsiov=""
  export crtmv=2.0.7
  export w3ev=_v2.1.0
  export w3nv=""
fi
#export gfsiov=${gfsiov:-_v1.1.0}
export crtmv=${crtmv:-2.0.7}
export crtmv_inc=${crtmv_inc:-v$crtmv}
export XMLPATH=${XMLPATH:-$NWPROD}
export BACIOPATH=${BACIOPATH:-$NWPROD/lib}
export xmlv=${xmlv:-""}
export w3ev=${w3ev:-_v2.0.3}
#export w3nv=${w3nv:-_v2.0.3}
export ipv=${ipv:-""}
export spv=${spv:-""}

if [ ${CLEAN:-YES}  = YES ] ; then make -f Makefile$BMPYXML clean ; fi

export make_nowrf=${make_nowrf:-NO}
if [ $make_post_lib = NO ] ; then
 if [ $make_post_exec = YES ] ; then
  if [ $make_nowrf = YES ] ; then
   make -f Makefile_nowrf${BMPYXML}
  else
   make -f Makefile${BMPYXML}
  fi
 fi
else
 if [ $make_post_exec = YES ] ; then
  if [ $make_nowrf = YES ] ; then
   make -f Makefile_nowrf${BMPYXML}
  else
   make -f Makefile${BMPYXML}
  fi
 fi
 export POSTLIBPATH=${POSTLIBPATH:-$(pwd)}
 if [ ${CLEAN:-YES}  = YES ] ; then rm -rf $POSTLIBPATH/incmod/post${BMPYXML}_4 ; fi
 mkdir -p $POSTLIBPATH/incmod/post${BMPYXML}_4
 make -f Makefile_lib${BMPYXML}
fi


