#!/bin/sh

#BSUB -o gdas_post.o%J
#BSUB -e gdas_post.o%J
#BSUB -J gdas_post
#BSUB -extsched 'CRAYLINUX[]'
#BSUB -W 01:00
#BSUB -q devhigh
#BSUB -P GFS-T2O
#BSUB -M 1000
#BSUB -cwd /gpfs/hps/emc/global/noscrub/emc.glopara/svn/gfs/work/gdas.v14.1.0/driver

set -x

export NODES=3
export ntasks=24
export ptile=8
export threads=1

# specify user's own post working directory for testing
export svndir=/gpfs/hps/emc/global/noscrub/emc.glopara/svn/gfs/work
export MP_LABELIO=yes

export OMP_NUM_THREADS=$threads


############################################
# Loading module
############################################
. $MODULESHOME/init/ksh
module load PrgEnv-intel ESMF-intel-haswell/3_1_0rp5 cfp-intel-sandybridge iobuf craype-hugepages2M craype-haswell
#module load cfp-intel-sandybridge/1.1.0
module use /gpfs/hps/nco/ops/nwprod/modulefiles
module load prod_envir
#module load prod_util
module load prod_util/1.0.4
module load grib_util/1.0.3

# specify PDY (the cycle start yyyymmdd) and cycle
export CDATE=2017052500
export PDY=`echo $CDATE | cut -c1-8`
export cyc=`echo $CDATE | cut -c9-10`
export cycle=t${cyc}z
export run_flag="hrly"   # if run_flag not "hrly", process "00 03 06 09"



# specify the directory environment for executable, it's either para or prod
export envir=prod

# set up running dir

export job=gdas_post_${cyc}
export pid=${pid:-$$}
export jobid=${job}.${pid}

export DATA=/gpfs/hps/stmp/$LOGNAME/test/$jobid
mkdir -p $DATA
cd $DATA
rm -f ${DATA}/*

####################################
# Specify RUN Name and model
####################################
export NET=gfs
#export RUN=gdas

####################################
# Determine Job Output Name on System
####################################
#export pgmout="OUTPUT.${pid}"
#export pgmerr=errfile

####################################
# SENDSMS  - Flag Events on SMS
# SENDCOM  - Copy Files From TMPDIR to $COMOUT
# SENDDBN  - Issue DBNet Client Calls
# RERUN    - Rerun posts from beginning (default no)
# VERBOSE  - Specify Verbose Output in global_postgp.sh
####################################
export SAVEGES=NO
export SENDSMS=NO
export SENDCOM=YES
export SENDDBN=NO
export RERUN=NO
export VERBOSE=YES

export HOMEglobal=${svndir}/global_shared.v14.1.0
export HOMEgfs=${svndir}/gfs.v14.1.0
export HOMEgdas=${svndir}/gdas.v14.1.0
                                                                                               
##############################################
# Define COM directories
##############################################
##export COMIN=$COMROOThps/gfs/para/gdas.${PDY}
export COMIN=/gpfs/hps/ptmp/emc.glopara/com2/gfs/para/gdas.${PDY}
# specify my own COMOUT dir to mimic operations
export COMOUT=/gpfs/hps/ptmp/$LOGNAME/com2/gfs/test/gdas.$PDY
mkdir -p $COMOUT

date

#export OUTTYP=4
# need to set FIXglobal to global share superstructure if testing post in non
# super structure environement
export FIXglobal=$svndir/global_shared.v14.1.0/fix
export APRUN="aprun -j 1 -n${ntasks} -N${ptile} -d${threads} -cc depth"
export nemsioget=$svndir/global_shared.v14.1.0/exec/nemsio_get

##export KEEPDATA=YES
export REMOVE_DATA=NO
#export POSTGRB2TBL=$HOMEglobal/parm/params_grib2_tbl_new
$HOMEgdas/jobs/JGDAS_NCEPPOST

#############################################################

date

echo $?



