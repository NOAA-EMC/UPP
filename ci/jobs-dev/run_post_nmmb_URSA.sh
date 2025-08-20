#!/bin/bash

#SBATCH -o out.post.nmmb
#SBATCH -e out.post.nmmb
#SBATCH -J nmmb_test
#SBATCH -t 00:20:00
#SBATCH -q batch
#SBATCH -N 7 --ntasks-per-node=4
#SBATCH -A ovp
#SBATCH --exclusive

set -x

# specify computation resources
export threads=1
export MP_LABELIO=yes
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
module list

msg="Starting nmmb test"
postmsg "$logfile" "$msg"

export POSTGPEXEC=${svndir}/exec/upp.x

# specify forecast start time and hour for running your post job
export startdate=2014120818
export fhr=03
export tmmark=tm00

# specify your running and output directory
export DATA=$rundir/nmmb_${startdate}
rm -rf $DATA; mkdir -p $DATA
cd $DATA

export NEWDATE=`$NDATE +${fhr} $startdate`
export YY=`echo $NEWDATE | cut -c1-4`
export MM=`echo $NEWDATE | cut -c5-6`
export DD=`echo $NEWDATE | cut -c7-8`
export HH=`echo $NEWDATE | cut -c9-10`

cat > itag <<EOF
&model_inputs
fileName='$homedir/data_in/nmmb/nmmb_hst_01_nio_00${fhr}h_00m_00.00s'
IOFORM='binarynemsio'
grib='grib2'
DateStr='${YY}-${MM}-${DD}_${HH}:00:00'
MODELNAME='NMM'
/
EOF

# copy fix data
cp $homedir/fix/nam_micro_lookup.dat ./eta_micro_lookup.dat
cp ${svndir}/parm/postxconfig-NT-NMM.txt ./postxconfig-NT.txt
cp ${svndir}/parm/params_grib2_tbl_new params_grib2_tbl_new

# Run the UPP
$APRUN ${POSTGPEXEC} < itag > outpost_nmmb_${NEWDATE}

mv BGDAWP${fhr}.${tmmark} BGDAWP${fhr}.${tmmark}.Grib2
mv BGRD3D${fhr}.${tmmark} BGRD3D${fhr}.${tmmark}.Grib2
mv BGRDSF${fhr}.${tmmark} BGRDSF${fhr}.${tmmark}.Grib2

################################################
# Compare with baseline data
################################################

# NMMB post processing generates 3 files
filelist="BGDAWP${fhr}.${tmmark}.Grib2 \
          BGRD3D${fhr}.${tmmark}.Grib2 \
          BGRDSF${fhr}.${tmmark}.Grib2"

for file in $filelist; do
export filein2=$file
ls -l ${filein2}
export err=$?

if [ $err = "0" ] ; then
 # use cmp to see if new pgb files are identical to the control one
 cmp ${filein2} $homedir/data_out_$compiler/nmmb/${filein2}.${machine}

 # if not bit-identical, use cmp_grib2_grib2 to compare each grib record
 export err1=$?
 if [ $err1 -eq 0 ] ; then
  msg="nmmb test: your new post executable generates bit-identical ${filein2} as the develop branch"
  echo $msg
 else
  msg="nmmb test: your new post executable did not generate bit-identical ${filein2} as the develop branch"
  echo $msg
  echo " start comparing each grib record and write the comparison result to *diff files"
  echo " check these *diff files to make sure your new post only change variables which you intend to change"
  $cmp_grib2_grib2 $homedir/data_out_$compiler/nmmb/${filein2}.${machine} ${filein2} > ${filein2}.diff
 fi
else
 msg="nmmb test: post failed using your new post executable to generate ${filein2}"
 echo $msg 2>&1 | tee -a TEST_ERROR
fi

postmsg "$logfile" "$msg"
done

echo "PROGRAM IS COMPLETE!!!!!" 2>&1 | tee SUCCESS
msg="Ending nmmb test"
postmsg "$logfile" "$msg"
