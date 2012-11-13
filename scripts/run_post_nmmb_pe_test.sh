#!/bin/sh

#
#@ account_no = NAM-MTN
#@ resources = ConsumableCpus(1) ConsumableMemory(1200)
#@ output = out.post_16pe.nmmb
#@ error = out.post_16pe.nmmb
#@ job_type = parallel
#@ class = dev
#@ class = 1
#@ group = class1onprod
#
#
#@ total_tasks = 16 
##@ task_affinity=cpu(1)
##@ blocking = unlimited
#@ node = 1
##@ node_usage = not_shared
#@ wall_clock_limit = 00:10:50
##@ preferences = Feature == "dev"
#@ network.MPI = csss,shared,us
#@ queue
#

set -x

# specify user's own post executable for testing
export POSTGPEXEC=/global/save/wx20hc/nceppost/sorc_nems/ncep_post
#export POSTGPEXEC=/meso/noscrub/wx20hc/post_RR_branch_new/ncep_post

# specify forecast start time and hour for running your post job
export startdate=2011062100
export fhr=12

# specify your running and output directory
export user=`whoami`
export DATA=/ptmp/${user}/post_nmmb_meso_${startdate}_16tasks

# specify your home directory 
export homedir=`pwd`/..

export tmmark=tm00
export MP_LABELIO=yes


mkdir -p $DATA
cd $DATA


#export allfhr="06"

#for fhr in $allfhr
#do

#export NEWDATE=$startdate
                                                                                       
export NEWDATE=`/nwprod/util/exec/ndate +${fhr} $startdate`
                                                                                       
export YY=`echo $NEWDATE | cut -c1-4`
export MM=`echo $NEWDATE | cut -c5-6`
export DD=`echo $NEWDATE | cut -c7-8`
export HH=`echo $NEWDATE | cut -c9-10`


cat > itag <<EOF
$homedir/data_in/nmm_b_history_nemsio.0${fhr}h_00m_00.00s
binarynemsio
${YY}-${MM}-${DD}_${HH}:00:00
NMM
EOF


rm -f fort.*

ln -sf $homedir/parm/nam_cntrl_cmaq.parm fort.14
ln -sf griddef.out fort.110
cp $homedir/parm/nam_micro_lookup.dat ./eta_micro_lookup.dat


${POSTGPEXEC} < itag > outpost_nems_${NEWDATE}


export err=$?
if [ $err = "0" ] ; then

 # operational NMMB post processing generates 3 files, start with BGDAWP first
 export filein=BGDAWP${fhr}.tm00
 # use cmp to see if new pgb files are identical to the control one
 cmp ${filein} $homedir/data_out/${filein}

 # if not bit-identical, use diffgb to compare each grib record
 export err1=$?
 if [ $err1 -eq 0 ] ; then
  echo " your new post executable generates bit-identical NMMB BGDAWP file as the trunk"
 else 
  echo " your new post executable did not generate bit-identical NMMB BGDAWP file as the trunk"
  echo " start comparing each grib record and write the comparison result to *diff files"
  echo " check these *diff files to make sure your new post only change variables which you intend to change"

  # since diffgb does not distinguish time stamp, seperate grib files by time stamp first
  grbtr ALL $filein ${filein}.'$TR'

  grbtr ALL $homedir/data_out/$filein ops.${filein}.'$TR'

  /global/save/wx20hc/bin/diffgb -x ops.${filein}.0 \
  ${filein}.0>${filein}.0.diff

  /global/save/wx20hc/bin/diffgb -x ops.${filein}.3 \
  ${filein}.3>${filein}.3.diff

  /global/save/wx20hc/bin/diffgb -x ops.${filein}.4 \
  ${filein}.4>${filein}.4.diff
 fi

 # next compare BGRD3D  
 export filein=BGRD3D${fhr}.tm00
 # use cmp to see if new pgb files are identical to the control one
 cmp ${filein} $homedir/data_out/${filein}

 # if not bit-identical, use diffgb to compare each grib record
 export err1=$?
 if [ $err1 -eq 0 ] ; then
  echo " your new post executable generates bit-identical NMMB BGRD3D file as the trunk"
 else 
  echo " your new post executable did not generate bit-identical NMMB BGRD3D file as the trunk"
  echo " start comparing each grib record and write the comparison result to *diff files"
  echo " check these *diff files to make sure your new post only change variables which you intend to change"

  # since diffgb does not distinguish time stamp, seperate grib files by time stamp first
  grbtr ALL $filein ${filein}.'$TR'

  grbtr ALL $homedir/data_out/$filein ops.${filein}.'$TR'

  /global/save/wx20hc/bin/diffgb -x ops.${filein}.0 \
  ${filein}.0>${filein}.0.diff

  /global/save/wx20hc/bin/diffgb -x ops.${filein}.3 \
  ${filein}.3>${filein}.3.diff

  /global/save/wx20hc/bin/diffgb -x ops.${filein}.4 \
  ${filein}.4>${filein}.4.diff
 fi
 
 # finally compare BGRDSF  
 export filein=BGRDSF${fhr}.tm00
 # use cmp to see if new pgb files are identical to the control one
 cmp ${filein} $homedir/data_out/${filein}

 # if not bit-identical, use diffgb to compare each grib record
 export err1=$?
 if [ $err1 -eq 0 ] ; then
  echo " your new post executable generates bit-identical NMMB BGRDSF file as the trunk"
 else 
  echo " your new post executable did not generate bit-identical NMMB BGRDSF file as the trunk"
  echo " start comparing each grib record and write the comparison result to *diff files"
  echo " check these *diff files to make sure your new post only change variables which you intend to change"

  # since diffgb does not distinguish time stamp, seperate grib files by time stamp first
  grbtr ALL $filein ${filein}.'$TR'

  grbtr ALL $homedir/data_out/$filein ops.${filein}.'$TR'

  /global/save/wx20hc/bin/diffgb -x ops.${filein}.0 \
  ${filein}.0>${filein}.0.diff

  /global/save/wx20hc/bin/diffgb -x ops.${filein}.3 \
  ${filein}.3>${filein}.3.diff

  /global/save/wx20hc/bin/diffgb -x ops.${filein}.4 \
  ${filein}.4>${filein}.4.diff
 fi
 
else

 echo "post failed using your new post executable"
 
fi  
  

#let "fhr=fhr+3"
#typeset -Z2 fhr

#done
