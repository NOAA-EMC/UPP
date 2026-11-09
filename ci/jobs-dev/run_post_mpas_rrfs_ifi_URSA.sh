#!/bin/bash 
 
#SBATCH -o out.post.mpas_rrfs_ifi
#SBATCH -e out.post.mpas_rrfs_ifi
#SBATCH -J mpas_rrfs_ifi_test
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

msg="Starting mpas_rrfs_ifi test"
postmsg "$logfile" "$msg"

export POSTGPEXEC=${svndir}/exec/upp_with_ifi.x

# specify forecast start time and hour for running your post job
export startdate=2026092806
export fhr=12

# specify your running and output directory
export DATA=$rundir/mpas_rrfs_ifi_${startdate}
rm -rf $DATA; mkdir -p $DATA
cd $DATA

export NEWDATE=`${NDATE} +${fhr} $startdate`
export YY=`echo ${NEWDATE} | cut -c1-4`
export MM=`echo ${NEWDATE} | cut -c5-6`
export DD=`echo ${NEWDATE} | cut -c7-8`
export HH=`echo ${NEWDATE} | cut -c9-10`

cat > itag <<EOF
&model_inputs
fileName='$homedir/data_in/mpas_rrfs_ifi/mpassit.${YY}-${MM}-${DD}_${HH}.00.00.nc'
fileNameFlux='$homedir/data_in/mpas_rrfs_ifi/mpassit.${YY}-${MM}-${DD}_${HH}.00.00.nc'
IOFORM='netcdfpara'
grib='grib2'
DateStr='${YY}-${MM}-${DD}_${HH}:00:00'
MODELNAME='RAPR'
SUBMODELNAME='MPAS'
/
&NAMPGB
numx=2
write_ifi_debug_files=.true.
/
EOF

# copy fix data
cp ${svndir}/parm/params_grib2_tbl_new params_grib2_tbl_new
cp ${svndir}/parm/postxconfig-NT-ifi.txt postxconfig-NT.txt
cp ${svndir}/fix/rap_micro_lookup.dat eta_micro_lookup.dat

# Run the UPP
${APRUN} ${POSTGPEXEC} < itag > outpost_mpas_rrfs_ifi_${NEWDATE}

################################################
# Compare with baseline data
################################################
fhr=`expr $fhr + 0`
fhr2=`printf "%02d" $fhr`

# MPAS RRFS_IFI post processing generates 1 file
filelist="IFIFIP.GrbF${fhr2}"

for file in $filelist; do
export filein2=$file
ls -l ${filein2}
export err=$?

if [ $err = "0" ] ; then
 # use cmp to see if new pgb files are identical to the control one
 cmp ${filein2} $homedir/data_out_$compiler/mpas_rrfs_ifi/${filein2}.${machine}

 # if not bit-identical, use cmp_grib2_grib2 to compare each grib record
 export err1=$?
 if [ $err1 -eq 0 ] ; then
  msg="mpas_rrfs_ifi test: your new post executable generates bit-identical ${filein2} as the develop branch"
  echo $msg
 else
  msg="mpas_rrfs_ifi test: your new post executable did not generate bit-identical ${filein2} as the develop branch"
  echo $msg
  echo " start comparing each grib record and write the comparison result to *diff files"
  echo " check these *diff files to make sure your new post only change variables which you intend to change"
  $cmp_grib2_grib2 $homedir/data_out_$compiler/mpas_rrfs_ifi/${filein2}.${machine} ${filein2} > ${filein2}.diff
 fi
else
 msg="mpas_rrfs_ifi test: post failed using your new post executable to generate ${filein2}"
 echo $msg 2>&1 | tee -a TEST_ERROR
fi

postmsg "$logfile" "$msg"
done

echo "PROGRAM IS COMPLETE!!!!!" 2>&1 | tee SUCCESS
msg="Ending mpas_rrfs_ifi test"
postmsg "$logfile" "$msg"
