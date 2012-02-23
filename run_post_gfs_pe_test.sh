#!/bin/sh

#
#@ account_no = GFS-MTN
#@ resources = ConsumableCpus(1) ConsumableMemory(1200)
#@ output = out.post_16pe.gfs
#@ error = out.post_16pe.gfs
#@ job_type = parallel
#@ class = dev
##@ node_usage = not_shared
#@ class = 1
#@ group = class1onprod
#
#
#@ total_tasks = 16 
##@ blocking = unlimited
#@ wall_clock_limit = 00:20:50
##@ preferences = Feature == "dev"
#@ network.MPI = csss,shared,us
#@ queue
#

set -x

# specify user's own post executable for testing
#export POSTGPEXEC=/global/save/wx20hc/nceppost/sorc_nems/ncep_post
export POSTGPEXEC=/global/save/wx20hc/RFC/post_upgrade_gfs_EnKF_2012/src/ncep_post
export POSTGPEXEC=/climate/save/wx20wa/esmf/nems/20120125_ngacgrb2/unipost/trunk/src/ncep_post

# specify PDY (the cycle start yyyymmdd) and cycle
export PDY=20110619
export cyc=00
export cycle=t${cyc}z

# specify your running and output directory
export user=`whoami`
export DATA=/ptmp/${user}/gfs.${PDY}_16pe

# this script mimics operational GFS post processing production
export MP_LABELIO=yes
# specify your home directory
export homedir=`pwd`/..


# specify the directory environment for executable, it's either para or prod
export envir=prod

mkdir -p $DATA
cd $DATA
#rm -f ${DATA}/*gfsio*

####################################
# Specify RUN Name and model
####################################
export NET=gfs
export RUN=gfs

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

####################################
# Specify Execution Areas
####################################
export GRBINDEX=/nwprod/util/exec/grbindex
                                                                                                         
export HOMEGLOBAL=/nw${envir}
export EXECGLOBAL=$HOMEGLOBAL/exec
export USHGLOBAL=$HOMEGLOBAL/ush
export FIXGLOBAL=$HOMEGLOBAL/fix
export PARMGLOBAL=$HOMEGLOBAL/parm
        
export HOMEUTIL=/nw${envir}/util
export EXECUTIL=$HOMEUTIL/exec
export FIXUTIL=$HOMEUTIL/fix
        
#export ERRSCRIPT=err_chk
#export LOGSCRIPT=startmsg
#export REDOUT='1>>'
#export REDERR='2>'

##############################################
# Define COM directories
##############################################
#export COMIN=/com/${NET}/${envir}/${RUN}.${PDY}
export COMIN=$homedir/data_in
#export COMOUT=/com/${NET}/${envir}/${RUN}.${PDY}
# specify my own COMOUT dir to mimic operations
export COMOUT=${DATA}
mkdir -p $COMOUT

##############################################
# Define GES directories
##############################################
gespath=/ptmp/${user}/GES
export GESdir=$gespath/${RUN}.${PDY}

####################################
# Specify Forecast Hour Range
####################################
export post_times=12
export SHOUR=192
export FHOUR=204
export FHINC=12
export HDMAX=36

####################################
# Specify Special Post Vars
####################################
export IGEN_ANL=81
export IGEN_FCST=96

export POSTGPVARS="KPO=47,PO=1000.,975.,950.,925.,900.,875.,850.,825.,800.,775.,750.,725.,700.,675.,650.,625.,600.,575.,550.,525.,500.,475.,450.,425.,400.,375.,350.,325.,300.,275.,250.,225.,200.,175.,150.,125.,100.,70.,50.,30.,20.,10.,7.,5.,3.,2.,1.,"
#######################################
# Specify Restart File Name to Key Off
#######################################
export restart_file=$COMIN/${RUN}.t${cyc}z.logf
                                                                                                         
####################################
# Specify Timeout Behavior of Post
#
# SLEEP_TIME - Amount of time to wait for
#              a restart file before exiting
# SLEEP_INT  - Amount of time to wait between
#              checking for restart files
####################################
export SLEEP_TIME=900
export SLEEP_INT=5

# use faster chgres
export CHGRESEXEC=$EXECGLOBAL/global_chgres
export CHGRESSH=$USHGLOBAL/global_chgres.sh
export CHGRESTHREAD=32

export ANOMGBEXEC=/nwprod/util/exec/anomgb
export FIXUTIL=/nwprod/util/fix
export COPYGB=/nwprod/util/exec/copygb
#############################################################
# specify new GFS control file that has more FD levels
export CTLFILE=$homedir/parm/gfs_cntrl.parm

# execute the script
export POSTGPSH=$homedir/scripts/global_nceppost.sh
$homedir/scripts/exgfs_nceppost.sh.sms.mytest

#############################################################


export err=$?
if [ $err = "0" ] ; then

 export filein=gfs.t${cyc}z.master.grbf${post_times}
 # use cmp to see if new pgb files are identical to the control one
 cmp ${filein} $homedir/data_out/${filein}

 # if not bit-identical, use diffgb to compare each grib record
 
 export err1=$?
 if [ $err1 -eq 0 ] ; then
  echo " your new post executable generates bit-identical GFS master file as the trunk"
 else 
  echo " your new post executable did not generate bit-identical master file as the trunk"
  echo " start comparing each grib record and write the comparison result to *diff files"
  echo " check these *diff files to make sure your new post only change variables which you intend to change"

  # since diffgb does not distinguish time stamp, seperate grib files by time stamp first
  grbtr ALL $filein ${filein}.'$TR'

  grbtr ALL $homedir/data_out/$filein ops.${filein}.'$TR'

  /global/save/wx20hc/bin/diffgb -x ops.${filein}.10 \
  ${filein}.10>${filein}.10.diff

  /global/save/wx20hc/bin/diffgb -x ops.${filein}.2 \
  ${filein}.2>${filein}.2.diff

  /global/save/wx20hc/bin/diffgb -x ops.${filein}.3 \
  ${filein}.3>${filein}.3.diff

  /global/save/wx20hc/bin/diffgb -x ops.${filein}.4 \
  ${filein}.4>${filein}.4.diff
 fi

 # compare simulated goes files as well

 export filein=gfs.t${cyc}z.special.grbf${post_times}
 # use cmp to see if new pgb files are identical to the control one
 cmp ${filein} $homedir/data_out/${filein}

 # if not bit-identical, use diffgb to compare each grib record

 export err1=$?
 if [ $err1 -eq 0 ] ; then
  echo " "
  echo " your new post executable generates bit-identical GFS special file as the trunk"
 else
  echo " your new post executable did not generate bit-identical special file as the trunk"
  echo " start comparing each grib record and write the comparison result to *diff files"
  echo " check these *diff files to make sure your new post only change variables which you intend to change"

  # since diffgb does not distinguish time stamp, seperate grib files
 by time stamp first
  grbtr ALL $filein ${filein}.'$TR'

  grbtr ALL $homedir/data_out/$filein ops.${filein}.'$TR'

  /global/save/wx20hc/bin/diffgb -x ops.${filein}.10 \
  ${filein}.10>${filein}.10.diff

  /global/save/wx20hc/bin/diffgb -x ops.${filein}.2 \
  ${filein}.2>${filein}.2.diff
  
  /global/save/wx20hc/bin/diffgb -x ops.${filein}.3 \
  ${filein}.3>${filein}.3.diff

  /global/save/wx20hc/bin/diffgb -x ops.${filein}.4 \
  ${filein}.4>${filein}.4.diff
 fi
fi 

echo $?
