#!/bin/ksh
set -x
mac=$(hostname | cut -c1-1)
mac2=$(hostname | cut -c1-2)
#export CLEAN=NO                                 # uncomment this if you don't want to clean before
                                                 # compiling
#
if [ $mac2 = ga ] ; then                         # For GAEA
 machine=gaea
 center=${center:-ncep}
elif [ $mac = z -o $mac = h -o $mac = f ] ; then # For ZEUS
 machine=zeus
elif [ $mac = t -o $mac = e -o $mac = g ] ; then # For WCOSS
 machine=wcoss
fi

if [ $machine = wcoss ] ; then
  export NETCDFPATH="/usrx/local/NetCDF/3.6.3"
  export WRFPATH="/nwprod/sorc/wrf_shared.v1.1.0"
  export NWPROD="/nwprod"
  export XMLPATH=$NWPROD
  export FC=mpiifort
  export CPP="/lib/cpp -P"
  export CPPFLAGS="-DLINUX"
# export OPTS="-O0 -openmp"
  export OPTS="-O3 -convert big_endian -fp-model source -openmp -xAVX"
  export LIST=""
  export FREE="-FR"
  export TRAPS=""
  export PROFILE=""
  export DEBUG=""
# export DEBUG="-g -check all -ftrapuv -convert big_endian -fp-stack-check -fstack-protector -heap-arrays -recursive -traceback"
elif [ $machine = zeus ] ; then
  export NETCDFPATH="/apps/netcdf/3.6.3/intel/lib"
  export WRFPATH="/scratch2/portfolios/NCEPDEV/meso/save/Dusan.Jovic/WRFV3"
  export NWPROD="/contrib/nceplibs/nwprod"
  export XMLPATH="/home/Hui-Ya.Chuang"
  export FC="ifort -lmpi -traceback"
  export CPP="/lib/cpp -P"
  export ARCH=""
  export CPPFLAGS="-DLINUX"
# export OPTS="-O0 -openmp -g"
  export export OPTS="-O3 -convert big_endian -traceback -g -fp-model source -openmp"
  export LIST=""
  export FREE="-FR"
  export TRAPS=""
  export PROFILE=""
  export DEBUG=""
# export DEBUG="-g -check all -ftrapuv -convert big_endian -fp-stack-check -fstack-protector -heap-arrays -recursive -traceback"
elif [ $machine = gaea ] ; then
  echo 'yet to be done'
  exit
  export NETCDFPATH="/apps/netcdf/3.6.3/intel/lib"
  export WRFPATH="/scratch2/portfolios/NCEPDEV/meso/save/Dusan.Jovic/WRFV3"
  export NWPROD="/contrib/nceplibs/nwprod"
  export XMLPATH="/home/Hui-Ya.Chuang"
  export FC="ifort -lmpi -traceback"
  export CPP="/lib/cpp -P"
  export ARCH=""
  export CPPFLAGS="-DLINUX"
# export OPTS=-O0-g
  export export OPTS="-O3 -convert big_endian -traceback -g -fp-model source"
  export LIST=""
  export FREE=-FR
  export TRAPS=""
  export PROFILE=""
  export DEBUG=""
# export DEBUG="-g -check all -ftrapuv -convert big_endian -fp-stack-check -fstack-protector -heap-arrays -recursive -traceback"
fi

if [ ${CLEAN:-YES}  = YES ] ; then make -f Makefile clean ; fi
make -f Makefile


