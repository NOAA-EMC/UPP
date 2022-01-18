#!/bin/bash
# Wen Meng 01/2020, Set up for cmake build.
# Wen Meng 01/2022, Add option for building with gtg code
############################################################

set -eu

usage() {
  echo
  echo "Usage: $0 [-g] [-w] -h"
  echo
  echo "  -p  installation prefix <prefix>    DEFAULT: ../install"
  echo "  -g  Build with GTG(users with gtg repos. access only)     DEFAULT: OFF"
  echo "  -w  Build without WRF-IO            DEFAULT: ON"
  echo "  -v  Build with cmake verbose        DEFAULT: NO"
  echo "  -h  display this message and quit"
  echo
  exit 1
}

prefix="../install"
gtg_opt=" -DBUILD_WITH_GTG=OFF"
wrfio_opt=" -DBUILD_WITH_WRFIO=ON"
verbose_opt=""
while getopts ":p:gwvh" opt; do
  case $opt in
    p)
      prefix=$OPTARG
      ;;
    g)
      gtg_opt=" -DBUILD_WITH_GTG=ON"
      ;;
    w)
      wrfio_opt=" -DBUILD_WITH_WRFIO=OFF"
      ;;
    v)
      verbose_opt="VERBOSE=1"
      ;;
    h|\?|:)
      usage
      ;;
  esac
done
cmake_opts=" -DCMAKE_INSTALL_PREFIX=$prefix"${wrfio_opt}${gtg_opt}

hostname
source ./detect_machine.sh
if [[ $(uname -s) == Darwin ]]; then
  readonly MYDIR=$(cd "$(dirname "$(greadlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)
else
  readonly MYDIR=$(cd "$(dirname "$(readlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)
fi
PATHTR=${PATHTR:-$( cd ${MYDIR}/.. && pwd )}

#Load required modulefiles
if [[ $MACHINE_ID != "unknown" ]]; then
   module purge
   module use $PATHTR/modulefiles
   modulefile=${MACHINE_ID}
   module load $modulefile
   module list
fi

rm -rf build install
mkdir build && cd build
cmake $cmake_opts ../..
make -j6 $verbose_opt 
make install

rm -rf $PATHTR/exec && mkdir $PATHTR/exec
cp $PATHTR/tests/install/bin/upp.x $PATHTR/exec/.
