#!/bin/bash 
 
#SBATCH -o out.post.hrrr_ifi
#SBATCH -e out.post.hrrr_ifi
#SBATCH -J hrrr_ifi_test
#SBATCH -t 00:30:00
#SBATCH -q batch
#SBATCH -A ovp
#SBATCH --exclusive
#SBATCH -N 2 --ntasks-per-node=24

# specify computation resources
export MP_LABELIO=yes
export threads=1
export OMP_NUM_THREADS=$threads
export APRUN="srun"

echo "starting time"
date

############################################
# Loading modules
############################################
module purge
module use $svndir/modulefiles
module load ursa_$compiler
module load wgrib2/3.6.0
module load prod_util/2.1.1
module load nccmp/1.9.1.0
module list

msg="Starting hrrr_ifi test"
postmsg "$logfile" "$msg"

export POSTGPEXEC=${svndir}/exec/upp_with_ifi.x

# specify forecast start time and hour for running your post job
export startdate=2025063004
export fhr=10

# specify your running and output directory
export DATA=$rundir/hrrr_ifi_${startdate}
rm -rf $DATA; mkdir -p $DATA
cd $DATA

export NEWDATE=`${NDATE} +${fhr} $startdate`
export YY=`echo ${NEWDATE} | cut -c1-4`
export MM=`echo ${NEWDATE} | cut -c5-6`
export DD=`echo ${NEWDATE} | cut -c7-8`
export HH=`echo ${NEWDATE} | cut -c9-10`

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
write_ifi_debug_files=.true.
/
EOF

# copy fix data
cp ${svndir}/parm/params_grib2_tbl_new params_grib2_tbl_new
cp ${svndir}/parm/postxconfig-NT-ifi.txt postxconfig-NT.txt
cp ${svndir}/fix/rap_micro_lookup.dat eta_micro_lookup.dat

# Run the UPP
${APRUN} ${POSTGPEXEC} < itag > outpost_hrrr_ifi_${NEWDATE}

################################################
# Compare with baseline data
################################################
fhr=`expr $fhr + 0`
fhr2=`printf "%02d" $fhr`

# HRRR_IFI post processing generates 1 file
filelist="IFIFIP.GrbF${fhr2}"

for file in $filelist; do
export filein2=$file
ls -l ${filein2}
export err=$?

if [ $err = "0" ] ; then
 # use cmp to see if new pgb files are identical to the control one
 cmp ${filein2} $homedir/data_out_$compiler/hrrr_ifi/${filein2}.${machine}

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
  $cmp_grib2_grib2 $homedir/data_out_$compiler/hrrr_ifi/${filein2}.${machine} ${filein2} > ${filein2}.diff
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
