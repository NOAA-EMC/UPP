#!/bin/bash
# Wen Meng 01/2020, Set up for cmake build.
# Wen Meng 01/2022, Add option for building with gtg code
# Sam Trahan 01/2023, Add option for building with libIFI
############################################################

set -eu

if [[ $(uname -s) == Darwin ]]; then
  readonly MYDIR=$(cd "$(dirname "$(greadlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)
else
  readonly MYDIR=$(cd "$(dirname "$(readlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)
fi
PATHTR=${PATHTR:-$( cd ${MYDIR}/.. && pwd )}
source ${PATHTR}/tests/detect_machine.sh

set_defaults() {
    delete_exec=YES
    upp_name="upp.x"
    load_ifi_module=NO
    prefix="../install"
    ifi_opt=" -DBUILD_WITH_IFI=OFF"
    build_ifi_executables_opt=" "
    build_ifi_executables=NO
    gtg_opt=" -DBUILD_WITH_GTG=OFF"
    nemsio_opt=" -DBUILD_WITH_NEMSIO=ON"
    wrfio_opt=" -DBUILD_WITH_WRFIO=ON"
    build_unit_tests_opt=" -DBUILD_TESTING=ON"
    more=" "
    verbose_opt=""
    debug_opt=""
    compiler="intel"
}

usage() {
  set_defaults # restore defaults so usage is correct
  echo
  echo "Usage: $0 [options]"
  echo
  echo "  -o  exe_name.x   Name of built UPP executable in exec. Default: $upp_name"
  echo "  -p  installation prefix <prefix>    DEFAULT: $prefix"
  echo "  -g  build with GTG(users with gtg repos. access only)     DEFAULT: ${gtg_opt#*=}"
  echo "  -i  build with libIFI(users with ifi install access only) DEFAULT: OFF"
  echo "  -I  build with libIFI (users with ifi repos. access only) DEFAULT: OFF"
  echo "  -B  build libIFI test programs (only valid with -I)       DEFAULT: OFF"
  echo "  -n  build with nemsio               DEFAULT: ${nemsio_opt#*=}"
  echo "  -w  build with WRF-IO               DEFAULT: ${wrfio_opt#*=}"
  echo "  -t  build unit tests                DEFAULT: ${build_unit_tests_opt#*=}"
  echo "  -v  build with cmake verbose        DEFAULT: OFF"
  echo "  -c  Compiler to use for build       DEFAULT: $compiler"
  echo "  -d  Debug mode of CMAKE_BUILD_TYPE  DEFAULT: Release"
  echo "  -a  Skip deletion of exec. Add new executables. DEFAULT: OFF"
  echo "  -Doption=value   Passes this option to cmake (can use more than once)"
  echo "  -h  display this message and quit"
  echo
  exit 1
}

set_defaults

while getopts ":p:gnwc:vhiIdBD:o:at" opt; do
  case $opt in
    a)
      delete_exec=NO
      ;;
    o)
      upp_name="$OPTARG"
      ;;
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
    t)
      build_unit_tests_opt=" -DBUILD_TESTING=OFF"
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

cmake_opts=" -DCMAKE_INSTALL_PREFIX=$prefix"${nemsio_opt}${wrfio_opt}${gtg_opt}${ifi_opt}${debug_opt}${build_ifi_executables_opt}${build_unit_tests_opt}${more}

#Load required modulefiles
if [[ $MACHINE_ID != "unknown" ]]; then
   if [ $MACHINE_ID == "wcoss2"  -o $MACHINE_ID == "wcoss2_a" ]; then
      module reset
   elif [ $MACHINE_ID == "container" ]; then
      source /usr/lmod/lmod/init/bash
      module purge
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
   modulefile=${MACHINE_ID}_${compiler}
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

# Provide host+compiler specific toolchains if available
CMAKE_TOOLCHAIN_FILE="${PATHTR}/cmake/toolchains/${MACHINE_ID}.${compiler}-toolchain.cmake"
if [[ -f "${CMAKE_TOOLCHAIN_FILE}" ]]; then
  cmake_opts="${cmake_opts} -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}"
fi

set -x
BUILD_DIR=${BUILD_DIR:-"build"}
rm -rf ${BUILD_DIR} install
mkdir -p ${BUILD_DIR} && cd ${BUILD_DIR}
cmake $cmake_opts ${PATHTR}
make -j${BUILD_JOBS:-6} $verbose_opt
make install

if [[ "$delete_exec" == YES ]] ; then
    rm -rf $PATHTR/exec
fi
test -d $PATHTR/exec || mkdir -p $PATHTR/exec
cp $prefix/bin/upp.x $PATHTR/exec/$upp_name
if [[ "$build_ifi_executables" == YES ]] ; then
    cp $prefix/bin/fip_runner $PATHTR/exec/.
fi
