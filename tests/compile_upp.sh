#!/bin/bash
# Wen Meng 01/2020, Set up for cmake build.
# Wen Meng 01/2022, Add option for building with gtg code
# Sam Trahan 01/2023, Add option for building with libIFI
############################################################

set -eu

usage() {
  echo
  echo "Usage: $0 [-p] [-g] [-w] [-v] [-c] [-i] [-d] -h"
  echo
  echo "  -p  installation prefix <prefix>    DEFAULT: ../install"
  echo "  -g  build with GTG(users with gtg repos. access only)     DEFAULT: OFF"
  echo "  -I  build with libIFI(users with ifi repos. access only)  DEFAULT: OFF"
  echo "  -i  build with libIFI(users with ifi install access only) DEFAULT: OFF"
  echo "  -w  build without WRF-IO            DEFAULT: ON"
  echo "  -v  build with cmake verbose        DEFAULT: NO"
  echo "  -c  Compiler to use for build       DEFAULT: intel"
  echo "  -d  Debug mode of CMAKE_BUILD_TYPE  DEFAULT: Release"
  echo "  -h  display this message and quit"
  echo
  exit 1
}

prefix="../install"
ifi_opt=" -DBUILD_WITH_IFI=OFF"
gtg_opt=" -DBUILD_WITH_GTG=OFF"
wrfio_opt=" -DBUILD_WITH_WRFIO=ON"
compiler="intel"
verbose_opt=""
debug_opt=""
while getopts ":p:gwc:vhiId" opt; do
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
    I)
      ifi_opt=" -DINTERNAL_IFI=ON"
      ;;
    i)
      ifi_opt=" -DREQUIRE_IFI=ON"
      ;;
    c)
      compiler=$OPTARG
      ;;
    v)
      verbose_opt="VERBOSE=1"
      ;;
    d)
      debug_opt=" -DCMAKE_BUILD_TYPE=Debug"
      ;;
    h|\?|:)
      usage
      ;;
  esac
done
cmake_opts=" -DCMAKE_INSTALL_PREFIX=$prefix"${wrfio_opt}${gtg_opt}${ifi_opt}${debug_opt}

if [[ $(uname -s) == Darwin ]]; then
  readonly MYDIR=$(cd "$(dirname "$(greadlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)
else
  readonly MYDIR=$(cd "$(dirname "$(readlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)
fi
PATHTR=${PATHTR:-$( cd ${MYDIR}/.. && pwd )}
source ${PATHTR}/tests/detect_machine.sh

#Load required modulefiles
if [[ $MACHINE_ID != "unknown" ]]; then
   if [ $MACHINE_ID == "wcoss2"  -o $MACHINE_ID == "wcoss2_a" ]; then
      module reset
   else
      module purge
   fi
   module use $PATHTR/modulefiles
   if [[ $compiler == "intel" ]]; then
      modulefile=${MACHINE_ID}
   else
      modulefile=${MACHINE_ID}_${compiler}
   fi
   if [ -f "${PATHTR}/modulefiles/${modulefile}" -o -f "${PATHTR}/modulefiles/${modulefile}.lua" ]; then
      echo "Building for machine ${MACHINE_ID}, compiler ${compiler}"
   else
      echo "Modulefile does not exist for machine ${MACHINE_ID}, compiler ${compiler}"
      exit 1
   fi
   module load $modulefile
   module list
fi

BUILD_DIR=${BUILD_DIR:-"build"}
rm -rf ${BUILD_DIR} install
mkdir -p ${BUILD_DIR} && cd ${BUILD_DIR}
cmake $cmake_opts ${PATHTR}
make -j${BUILD_JOBS:-6} $verbose_opt
make install

rm -rf $PATHTR/exec && mkdir -p $PATHTR/exec
cp $prefix/bin/upp.x $PATHTR/exec/.
