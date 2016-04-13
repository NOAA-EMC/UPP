SHELL=/bin/sh

####################################################################################################
#
# post using module compile standard
#
# 10/15 Lin Gan:        Create module load version
# 01/16 Lin Gan:	Update to use GFS Vertical Structure
#
#####################################################################################################
#####################################################################################################


# Lin Gan Module Load
module purge

# Lin Gan modifiy to use NCO vertical structure prefix for NCO deployment - 20160131
module load ../modulefiles/post/v7.0.0
module list

cd ncep_post.fd
make -f makefile_wcoss_module clean
make -f makefile_wcoss_module 

cp ncep_post ../../exec/
