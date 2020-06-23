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
validmachines=(theia jet wcoss_dell_p3 wcoss cray-intel hera orion odin stampede)

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
   if [ $mac2 = tf ] ; then                         # For Theia
      machine=theia
   elif [ $mac = f  ] ; then                        # For Jet 
      machine=jet
   elif [ $mac = v -o $mac = m  ] ; then            # For Dell
      machine=wcoss_dell_p3
   elif [ $mac = t -o $mac = e -o $mac = g ] ; then # For WCOSS
      machine=wcoss
   elif [ $mac = l -o $mac = s ] ; then             #    wcoss_c (i.e. luna and surge)
      export machine=cray-intel
   elif [ $mac2 = hf ] ; then                       # For Hera
      machine=hera
   elif [ $mac = O ] ; then
      machine=orion
   elif [ $mac2 = od ] ; then
      machine=odin
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
orion)                                 # For Orion
 . /etc/profile
 ;;
odin)                                  # For Odin at NSSL
 . /etc/profile
 . /etc/profile.d/modules.sh
 ;;
stampede)
 module purge
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

exit 0
