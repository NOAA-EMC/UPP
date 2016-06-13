SHELL=/bin/sh

####################################################################################################
#
# post using module compile standard
#
# 10/15 Lin Gan:        Create module load version
#
#####################################################################################################
#####################################################################################################


# Lin Gan Module Load
module purge

module load ../../modulefiles/post/v7.0.0

make -f makefile_wcoss_module clean
make -f makefile_wcoss_module 
