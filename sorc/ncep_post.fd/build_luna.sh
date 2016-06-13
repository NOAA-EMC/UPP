#SHELL=/bin/sh

####################################################################################################
#
# post using module compile standard
#
# 10/15 Lin Gan:        Create module load version
# 03/16 Lin Gan:	Fix module issue with NCO, compiler option for bit identical results
#
#####################################################################################################
#####################################################################################################


# Lin Gan Module Load
module purge

module load -a ../../modulefiles/post/v7.0.0-cray-intel

module list

make -f makefile_luna_module clean
make -f makefile_luna_module
export POSTLIBPATH=${POSTLIBPATH:-$(pwd)}
mkdir -p include/post_4
make -f Makefile_lib
cp ncep_post ../../exec
