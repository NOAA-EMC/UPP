#!/bin/ksh
set -x
mac=$(hostname | cut -c1-1)
mac2=$(hostname | cut -c1-2)
 export CLEAN=NO                                 # uncomment this if you don't want to clean
                                                 # before compiling
#debug=YES                                       # to turn on debug mode
 make_post_lib=YES                               # to create post library
#
if [ $mac2 = ga ] ; then                         # For GAEA
 machine=gaea
 center=${center:-ncep}
elif [ $mac = z -o $mac = h -o $mac = f ] ; then # For ZEUS
 machine=zeus
elif [ $mac = t -o $mac = e -o $mac = g ] ; then # For WCOSS
 machine=wcoss
fi
debug=${debug:-NO}
make_post_lib=${make_post_lib:-NO}
if [ $machine = wcoss ] ; then
  export NETCDFPATH="/usrx/local/NetCDF/3.6.3"
  export WRFPATH="/nwprod/sorc/wrf_shared.v1.1.0"
  export NWPROD="/nwprod"
  export XMLPATH=$NWPROD
  export IPPATH=$NWPROD
  export SPPATH=/usrx/local/nceplibs
  export ipv=""
  export spv=_v2.0.2p
  export crtmv=2.0.6
  export crtmv_inc=$crtmv
  export xmlv=_v2.0.0
  export FC=mpiifort
  export CPP="/lib/cpp -P"
  export CPPFLAGS="-DLINUX"
  export CC=cc
  if [ $debug = YES ] ; then
    export OPTS="-O0 -openmp "
    export DEBUG="-g -check all -ftrapuv -convert big_endian -fp-stack-check -fstack-protector -heap-arrays -recursive -traceback"
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
  export FC="ifort -lmpi -traceback"
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
# export FC="ifort -lmpi -traceback"
  export FC="ftn -traceback"
  export CPP="/lib/cpp -P"
  export ARCH=""
  export CPPFLAGS="-DLINUX"
  export CC=icc
  if [ $debug = YES ] ; then
    export OPTS=-O0-g
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
export xmlv=${xmlv:-""}
export w3ev=${w3ev:-_v2.0.3}
#export w3nv=${w3nv:-_v2.0.3}
export ipv=${ipv:-""}
export spv=${spv:-""}

if [ ${CLEAN:-YES}  = YES ] ; then make -f Makefile clean ; fi

if [ $make_post_lib = NO ] ; then
 make -f Makefile
else
 export POSTLIBPATH=${POSTLIBPATH:-$(pwd)}
 if [ ${CLEAN:-YES}  = YES ] ; then rm -rf $POSTLIBPATH/incmod/post_4 ; fi
 mkdir -p $POSTLIBPATH/incmod/post_4
 make -f Makefile_lib
fi


