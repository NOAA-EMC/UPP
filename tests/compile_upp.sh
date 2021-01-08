#!/bin/bash
# Wen Meng 01/2020, Set up for cmake build.
#############################################

#set -eu

#List of valid machines:
validmachines=(wcoss_dell_p3 cray-intel hera orion jet)

function usage {
   echo "Usage:"
   echo " $0 machinename"
   echo ""
   echo " Valid values for 'machinename' are: ${validmachines[@]}"
   exit 1
}

if [ "$#" -eq 0 ]; then

  #Check to see if we are building the old way
  module purge

  set -x
  mac=$(hostname | cut -c1-1)
  mac2=$(hostname | cut -c1-2)

  if [ $mac2 = hf ] ; then   # For Hera
    machine=hera
  elif [ $mac = v -o $mac = m  ] ; then  # For WCOSS Dell
    machine=wcoss_dell_p3
  elif [ $mac = O ] ; then  # For orion
    machine=orion
  elif [ $mac = l -o $mac = s ] ; then  # For WCOSS Cray
    machine=cray-intel
  else
      echo ""
      echo "ERROR ERROR ERROR"
      echo ""
      echo "Error: To use this build script without arguments you must be on a valid machine"
      echo "Valid machines are:"
      echo "${validmachines[@]}"
      echo ""
      echo "ERROR ERROR ERROR"
  fi
elif [ "$#" -gt 1 ]; then
  echo "Error: too many input arguments"
  exit 2
else
  machine=$1 
fi

case $machine in
jet)                                   # For Jet
 module purge
 . /etc/profile
 . /etc/profile.d/modules.sh
 ;;
wcoss_dell_p3)                         # For Dell
 module purge
 . $MODULESHOME/init/bash
 ;;
cray-intel)                            # For wcoss_c (i.e. luna and surge)
 module purge
 ;;
hera)                                  # For Hera
 . /etc/profile
 . /etc/profile.d/modules.sh
 ;;
orion)                                 # For Orion
 . /etc/profile
 ;;
*)
 set +x
 echo "ERROR: Invalid machine name specified"
 usage
 ;;
esac

#Load requited modules
moduledir=`dirname $(readlink -f ../modulefiles)`
echo $moduledir
module use ${moduledir}/modulefiles
module load ${machine}
module list

set -x

rm -rf build install
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=../install -DBUILD_WITH_WRFIO=ON ../..
make -j12
make install

cd ../install
rm -rf ../../exec && mkdir ../../exec
cp bin/upp.x ../../exec/.
