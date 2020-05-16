SHELL=/bin/sh

####################################################################################################
#
# post using module compile standard
#
# 10/15 Lin Gan:        Create module load version
# 01/16 Lin Gan:	Update to use GFS Vertical Structure
# 07/16 J. Carley:      Generalize for other machines using modules
# 07/18 Wen Meng:       Set post to v8.0.0 for fv3gfs
#
#####################################################################################################
#####################################################################################################


# Lin Gan Module Load
module purge
set -x
mac=$(hostname | cut -c1-1)
mac2=$(hostname | cut -c1-2)
if [ $mac2 = tf ] ; then                        # For Theia
 machine=theia
 . /etc/profile
 . /etc/profile.d/modules.sh
elif [ $mac = f  ] ; then                       # For Jet
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
elif [ $mac2 = hf ] ; then                        # For Hera
 machine=hera
 . /etc/profile
 . /etc/profile.d/modules.sh
elif [ $mac = O ] ; then           # For Orion
 machine=orion
 . /etc/profile
elif [ $mac2 = od ] ; then                        # For Odin at NSSL
 machine=odin
 . /etc/profile
 . /etc/profile.d/modules.sh
fi

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
