SHELL=/bin/sh

module purge
set -x
mac=$(hostname | cut -c1-1)
mac2=$(hostname | cut -c1-2)

if [ $mac2 = hf ] ; then                        # For Hera
 machine=hera
 . /etc/profile
 . /etc/profile.d/modules.sh
elif [ $mac = f  ] ; then            # For Jet
 machine=jet
 . /etc/profile
 . /etc/profile.d/modules.sh
elif [ $mac = v -o $mac = m  ] ; then            # For Dell
 machine=wcoss_dell_p3
 . $MODULESHOME/init/bash
elif [ $mac = t -o $mac = e -o $mac = g ] ; then # For WCOSS
 machine=wcoss
 . /usrx/local/Modules/default/init/bash
elif [ $mac = l -o $mac = s ] ; then             #    wcoss_c (i.e. luna and surge)
 export machine=cray-intel
elif [ $mac = O ] ; then           # For Orion
 machine=orion
 . /etc/profile
fi
export version=${1:-"v8.0.0"}

moduledir=`dirname $(readlink -f ../../modulefiles/post)`
module use -a ${moduledir}
module load post/lib-${machine}
#module load nceppost_modulefile

#
module list

#sleep 1

BASE=`pwd`

#####################################
cd ${BASE}
rm *.o *.mod  incmod
#mkdir -m 775 -p $BASE/../../lib/include/ncep_post_${version}_4
make -f makefile_lib clean
mkdir -m 775 -p include/ncep_post_4
make -f makefile_lib

exit 0

