#!/bin/bash
# Wen Meng 01/2020, Set up for cmake build.
#############################################

set -x

#Clean loaded modules
module purge

hostname
source ./detect_machine.sh
if [[ $(uname -s) == Darwin ]]; then
  readonly MYDIR=$(cd "$(dirname "$(greadlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)
else
  readonly MYDIR=$(cd "$(dirname "$(readlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)
fi
PATHTR=${PATHTR:-$( cd ${MYDIR}/.. && pwd )}

#Load required modulefiles
module use $PATHTR/modulefiles
modulefile=${MACHINE_ID}
module load $modulefile
module list

rm -rf build install
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=../install -DBUILD_WITH_WRFIO=ON ../..
make -j6 
make install

rm -rf $PATHTR/exec && mkdir $PATHTR/exec
cp $PATHTR/tests/install/bin/upp.x $PATHTR/exec/.
