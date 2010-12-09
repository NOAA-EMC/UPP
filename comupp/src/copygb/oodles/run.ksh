#!/bin/ksh
#####################################################
# machine set up (users should change this part)
#####################################################
#
#$ -S /bin/ksh
#$ -N mycopygb
#$ -cwd
#$ -r y
## #$ -pe ncomp 1
#$ -pe serial 1
#$ -l h_rt=0:10:00,h_vmem=9.0G 
## #$ -A wrfruc
#$ -A rtrr

#
# GSIPROC = processor number used for GSI analysis
#------------------------------------------------
ulimit -s 512000

copygb -g"255 3 451 337 16281 -126138 8 -95000 13545 13545 0 64 25000 25000" -N copygb_nat.parm -x wrfnat00 wrfnat_130_00.tm00
## copygb -g"255 3 451 337 16281 -126138 8 -95000 13545 13545 0 64 25000 25000" -N copygb_nat.parm -x wrfnat_rr_00.grib1 wrfnat_130_00.tm00


exit
