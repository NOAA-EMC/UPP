#!/bin/bash

####################################################################################################
#
# post using module compile standard
#
# 10/15 Lin Gan:        Create module load version
# 01/16 Lin Gan:	Update to use GFS Vertical Structure
# 07/16 J. Carley:      Generalize for other machines using modules
# 07/18 Wen Meng:       Set post to v8.0.0 for fv3gfs
# 10/19 M Kavulich:     Provide machine name as an input argument
#
#####################################################################################################
#####################################################################################################

#List of valid machines:
validmachines=(theia jet wcoss_dell_p3 wcoss cray-intel hera)

function usage {
   echo "Usage:"
   echo " $0 machinename"
   echo ""
   echo " Valid values for 'machinename' are: ${validmachines[@]}"
   exit 1
}

if [ "$#" -eq 0 ]; then
   echo "Error: To use this build script you must provide the machine name"
   usage
elif [ "$#" -gt 1 ]; then
   echo "Error: too many input arguments"
   exit 1
fi

machine=$1

# Lin Gan Module Load
set -x
case $machine in
theia)                                 # For Theia
 module purge
 . /etc/profile
 . /etc/profile.d/modules.sh
 ;;
jet)                                   # For Jet
 module purge
 . /etc/profile
 . /etc/profile.d/modules.sh
 ;;
wcoss_dell_p3)                         # For Dell
 module purge
 . $MODULESHOME/init/bash                 
 ;;
wcoss)                                 # For WCOSS
 module purge
 . /usrx/local/Modules/default/init/bash
 ;;
cray-intel)                            # For wcoss_c (i.e. luna and surge)
 module purge
 ;;
hera)                                  # For Hera
 . /etc/profile
 . /etc/profile.d/modules.sh
 ;;
*)
 set +x
 echo "ERROR: Invalid machine name specified"
 usage
 ;;
esac

# Lin Gan modifiy to use NCO vertical structure prefix for NCO deployment - 20160131
moduledir=`dirname $(readlink -f ../modulefiles/post)`
module use ${moduledir}
module load post/v8.0.0-${machine}
module list

cd ncep_post.fd

make -f makefile_module clean
make -f makefile_module

if [ ! -d "../../exec" ] ; then
  mkdir -p ../../exec
fi
cp ncep_post ../../exec/
