#!/bin/sh
set -x
export FCMP=${1:-ifort}
export LIBDIR=${2:-${LIBDIR:-/nwprod/lib}}
mac=$(hostname | cut -c1-1)

if [ $mac = f ] ; then  #For Zeus
  export LIBDIR=/contrib/nceplibs/nwprod/lib
  export LIBSMOD="-L$LIBDIR -lw3nco_4 -lbacio_4"
else
  export LIBSMOD="-L$LIBDIR -lw3nco_4 -lbacio_4"
fi
export FFLAGSM=""
export FFLAGS=""
make -f makefile

