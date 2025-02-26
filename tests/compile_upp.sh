#!/bin/bash
# Wen Meng 01/2020, Set up for cmake build.
# Wen Meng 01/2022, Add option for building with gtg code
# Sam Trahan 01/2023, Add option for building with libIFI
############################################################

set -eu

usage() {
  echo
  echo "Usage: $0 [options]"
  echo
  echo "  -p  installation prefix <prefix>    DEFAULT: ../install"
  echo "  -g  build with GTG(users with gtg repos. access only)     DEFAULT: OFF"
  echo "  -i  build with libIFI(users with ifi install access only) DEFAULT: OFF"
  echo "  -I  build with libIFI (users with ifi repos. access only) DEFAULT: OFF"
  echo "  -B  build libIFI test programs (only valid with -I)       DEFAULT: OFF"
  echo "  -n  build without nemsio            DEFAULT: ON"
  echo "  -w  build without WRF-IO            DEFAULT: ON"
  echo "  -v  build with cmake verbose        DEFAULT: NO"
  echo "  -c  Compiler to use for build       DEFAULT: intel"
  echo "  -d  Debug mode of CMAKE_BUILD_TYPE  DEFAULT: Release"
  echo "  -Doption=value   Passes this option to cmake (can use more than once)"
  echo "  -h  display this message and quit"
  echo
  exit 1
}

load_ifi_module=NO
prefix="../install"
ifi_opt=" -DBUILD_WITH_IFI=OFF"
build_ifi_executables_opt=" "
build_ifi_executables=NO
gtg_opt=" -DBUILD_WITH_GTG=OFF"
nemsio_opt=" -DBUILD_WITH_NEMSIO=ON"
wrfio_opt=" -DBUILD_WITH_WRFIO=ON"
more=" "
compiler="intel"
verbose_opt=""
debug_opt=""
while getopts ":p:gnwc:vhiIdBD:" opt; do
  case $opt in
    D)
      more="$more -$opt$OPTARG"
      ;;
    p)
      prefix=$OPTARG
      ;;
    g)
      gtg_opt=" -DBUILD_WITH_GTG=ON"
      ;;
    B)
      build_ifi_executables_opt=" -DBUILD_IFI_EXECUTABLES=ON"
      build_ifi_executables=YES
      ;;
    n)
      nemsio_opt=" -DBUILD_WITH_NEMSIO=OFF"
      ;;
    w)
      wrfio_opt=" -DBUILD_WITH_WRFIO=OFF"
      ;;
    I)
      ifi_opt=" -DINTERNAL_IFI=ON"
      ;;
    i)
      ifi_opt=" -DREQUIRE_IFI=ON"
      load_ifi_module=YES
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

if [[ ! -z $debug_opt && $ifi_opt =~ INTERNAL.*=ON ]] ; then
    echo ENABLING IFI DEBUG
    # When building debug mode with internal IFI, also enable debugging in IFI.
    # This includes bounds checking in much of the libIFI C++ library.
    debug_opt="$debug_opt -DIFI_DEBUG=ON"
fi

cmake_opts=" -DCMAKE_INSTALL_PREFIX=$prefix"${nemsio_opt}${wrfio_opt}${gtg_opt}${ifi_opt}${debug_opt}${build_ifi_executables_opt}${more}

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
   elif [[ "$MACHINE_ID" =~ gaea* ]] ; then
       module reset
       # Unset the read-only variables $PELOCAL_PRGENV and $RCLOCAL_PRGENV
       gdb -ex 'call (int) unbind_variable("PELOCAL_PRGENV")' \
           -ex 'call (int) unbind_variable("RCLOCAL_PRGENV")' \
           --pid=$$ --batch
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
   if [[ "$load_ifi_module" == YES ]] ; then
       echo "Loading modulefile for external libIFI library"
       module load ${modulefile}_external_ifi
   fi
   if [[ "$build_ifi_executables" == YES ]] ; then
       echo "Loading libIFI executables' prerequisites"
       module load ${modulefile}_ifi_test_prereqs
   fi
   module list
fi

set -x
BUILD_DIR=${BUILD_DIR:-"build"}
rm -rf ${BUILD_DIR} install
mkdir -p ${BUILD_DIR} && cd ${BUILD_DIR}
cmake $cmake_opts ${PATHTR}
make -j${BUILD_JOBS:-6} $verbose_opt
make install

rm -rf $PATHTR/exec && mkdir -p $PATHTR/exec
cp $prefix/bin/upp.x $PATHTR/exec/.
if [[ "$build_ifi_executables" == YES ]] ; then
    cp $prefix/bin/fip2-lookalike.x $PATHTR/exec/.
fi
