#!/bin/bash 
 
#PBS -o out.post.hrrr_ifi
#PBS -e out.post.hrrr_ifi
#PBS -N hrrr_ifi_test
#PBS -l walltime=00:30:00
#PBS -q debug
#PBS -A GFS-DEV
#PBS -l place=vscatter,select=2:ncpus=24
#PBS -V

set -x

# specify computation resource
export threads=1
export OMP_NUM_THREADS=$threads
export APRUN="mpiexec -l -n 48 -ppn 24"

echo "starting time"
date

######################################################################
# Purpose: to run RAP post processing
######################################################################
module reset
module load intel/19.1.3.304
module load PrgEnv-intel/8.1.0
module load craype/2.7.8
module load cray-mpich/8.1.7
module load cray-pals/1.0.12
module load hdf5/1.10.6
module load netcdf/4.7.4
module load libjpeg/9c
module load prod_util/2.0.8
module load crtm/2.4.0.1
module load wgrib2/2.0.8
module list

msg="Starting hrrr_ifi test"
postmsg "$logfile" "$msg"


export POSTGPEXEC=${svndir}/exec/upp.x

# CALL executable job script here

# specify your running and output directory
export startdate=2020060118
export fhr=04
export DATA=$rundir/hrrr_ifi_${startdate}

export NEWDATE=`${NDATE} +${fhr} $startdate`

export YY=`echo ${NEWDATE} | cut -c1-4`
export MM=`echo ${NEWDATE} | cut -c5-6`
export DD=`echo ${NEWDATE} | cut -c7-8`
export HH=`echo ${NEWDATE} | cut -c9-10`

rm -rf $DATA; mkdir -p $DATA
cd $DATA

cat > itag <<EOF
&model_inputs
fileName='$homedir/data_in/hrrr/wrfout_d01_${YY}-${MM}-${DD}_${HH}_00_00'
IOFORM='netcdf'
grib='grib2'
DateStr='${YY}-${MM}-${DD}_${HH}:00:00'
MODELNAME='RAPR'
/
&NAMPGB
KPO=47,PO=2.,5.,7.,10.,20.,30.,50.,70.,75.,100.,125.,150.,175.,200.,225.,250.,275.,300.,325.,350.,375.,400.,425.,450.,475.,500.,525.,550.,575.,600.,625.,650.,675.,700.,725.,750.,775.,800.,825.,850.,875.,900.,925.,950.,975.,1000.,1013.2
write_ifi_debug_files=.false.
/
EOF
#FMIN

#copy xml
cp ${svndir}/parm/params_grib2_tbl_new params_grib2_tbl_new
cp ${svndir}/parm/postxconfig-NT-ifi.txt postxconfig-NT.txt
cp ${svndir}/fix/rap_micro_lookup.dat eta_micro_lookup.dat

${APRUN} ${POSTGPEXEC} < itag > wrfpost2.out

# operational hrrr post processing generates 3 files
filelist="IFIFIP.GrbF04"

for file in $filelist; do
export filein2=$file
ls -l ${filein2}
export err=$?

if [ $err = "0" ] ; then

 # operational hrrr post processing generates 3 files, start with BGDAWP first
 # use cmp to see if new pgb files are identical to the control one
 cmp ${filein2} $homedir/data_out/hrrr_ifi/${filein2}.${machine}

 # if not bit-identical, use cmp_grib2_grib2 to compare each grib record
 export err1=$?
 if [ $err1 -eq 0 ] ; then
  msg="hrrr_ifi test: your new post executable generates bit-identical ${filein2} as the develop branch"
  echo $msg
 else
  msg="hrrr_ifi test: your new post executable did not generate bit-identical ${filein2} as the develop branch"
  echo $msg
  echo " start comparing each grib record and write the comparison result to *diff files"
  echo " check these *diff files to make sure your new post only change variables which you intend to change"
  $cmp_grib2_grib2 $homedir/data_out/hrrr_ifi/${filein2}.${machine} ${filein2} > ${filein2}.diff
 fi


else

 msg="hrrr_ifi test: post failed using your new post executable to generate ${filein2}"
 echo $msg 2>&1 | tee -a TEST_ERROR

fi
postmsg "$logfile" "$msg"
done

echo "PROGRAM IS COMPLETE!!!!!" 2>&1 | tee SUCCESS
msg="Ending hrrr_ifi test"
postmsg "$logfile" "$msg"
