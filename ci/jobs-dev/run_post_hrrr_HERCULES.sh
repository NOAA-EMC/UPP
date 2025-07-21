#!/bin/bash 
 
#SBATCH -o out.post.hrrr
#SBATCH -e out.post.hrrr
#SBATCH -J hrrr_test
#SBATCH -t 00:20:00
#SBATCH -q batch
#SBATCH -A nems
#SBATCH -N 2 --ntasks-per-node=24

set -x

# specify computation resources
export MP_LABELIO=yes
export threads=3
export OMP_NUM_THREADS=$threads
export APRUN="srun"

echo "starting time"
date

############################################
# Loading modules
############################################
module use ${svndir}/modulefiles
module load hercules_$compiler
module load prod_util/2.1.1
module load wgrib2/3.6.0
module list

ulimit -s unlimited

msg="Starting hrrr test"
postmsg "$logfile" "$msg"

export POSTGPEXEC=${svndir}/exec/upp.x     

# specify forecast start time and hour
export startdate=2025063004
export fhr=10

# specify your running and output directory
export DATA=$rundir/hrrr_${startdate}
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
/
EOF

# copy fix data
cp $homedir/fix/raphrrr_fix/* .
cp ${svndir}/parm/params_grib2_tbl_new params_grib2_tbl_new
cp ${svndir}/parm/postxconfig-NT-hrrr.txt postxconfig-NT.txt

# Run the UPP
${APRUN} ${POSTGPEXEC} < itag > outpost_hrrr_${NEWDATE}

################################################
# Compare with baseline data
################################################
fhr=`expr $fhr + 0`
fhr2=`printf "%02d" $fhr`

# HRRR post processing generates 3 files
filelist="WRFTWO.GrbF${fhr2} \
          WRFPRS.GrbF${fhr2} \
          WRFNAT.GrbF${fhr2}"

for file in $filelist; do
export filein2=$file
ls -l ${filein2}
export err=$?

if [ $err = "0" ] ; then
 # use cmp to see if new pgb files are identical to the control one
 cmp ${filein2} $homedir/data_out_$compiler/hrrr/${filein2}.${machine}

 # if not bit-identical, use cmp_grib2_grib2 to compare each grib record
 export err1=$?
 if [ $err1 -eq 0 ] ; then
  msg="hrrr test: your new post executable generates bit-identical ${filein2} as the develop branch"
  echo $msg
 else
  msg="hrrr test: your new post executable did not generate bit-identical ${filein2} as the develop branch"
  echo $msg
  echo " start comparing each grib record and write the comparison result to *diff files"
  echo " check these *diff files to make sure your new post only change variables which you intend to change"
  $cmp_grib2_grib2 $homedir/data_out_$compiler/hrrr/${filein2}.${machine} ${filein2} > ${filein2}.diff
 fi
else
 msg="hrrr test: post failed using your new post executable to generate ${filein2}"
 echo $msg 2>&1 | tee -a TEST_ERROR
fi

postmsg "$logfile" "$msg"
done

echo "PROGRAM IS COMPLETE!!!!!" 2>&1 | tee SUCCESS
msg="Ending hrrr test"
postmsg "$logfile" "$msg"
