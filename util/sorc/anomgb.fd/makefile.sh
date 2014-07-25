#!/bin/sh
set -x
export FCMP=${1:-ifort}
export LIBDIR=${2:-${LIBDIR:-/nwprod/lib}}
mac=$(hostname | cut -c1-1)

if [ $mac = f ] ; then  #For Zeus
  export LIBDIR=/contrib/nceplibs/nwprod/lib
  export LIBSMOD="-L$LIBDIR -lw3nco_4 -lip_d -lsp_d -lbacio_4"
else
  export LIBSMOD="-L$LIBDIR -lw3nco_4 -lip_d -lsp_d -lbacio_4"
fi
export FFLAGSM="-O3 -g -convert big_endian -assume byterecl -assume noold_ldout_format -r8"
export FFLAGS="-mkl"
make -f makefile
