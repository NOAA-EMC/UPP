#!/bin/bash

#PBS -o out.post.nmmb
#PBS -e out.post.nmmb
#PBS -N nmmb.test
#PBS -l walltime=00:30:00
#PBS -q debug
#PBS -A GFS-DEV
#PBS -l place=vscatter,select=2:ncpus=12
#PBS -V

set -x

# specify computation resource
export threads=1
export OMP_NUM_THREADS=$threads
export APRUN="mpiexec -l -n 24 -ppn 12" 

echo "starting time" 
date

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
module load wgrib2/2.0.8
module list


msg="Starting nmmb test"
postmsg "$logfile" "$msg"


# specify user's own post executable for testing
export POSTGPEXEC=${svndir}/exec/upp.x


# specify forecast start time and hour for running your post job
export startdate=2014120818
export fhr=03

# specify your running and output directory
export DATA=$rundir/nmmb_${startdate}

# specify your home directory 
#export homedir=`pwd`/..

export tmmark=tm00

rm -rf $DATA; mkdir -p $DATA
cd $DATA

echo $homedir
echo $NDATE
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

rm -f fort.*

cp $homedir/fix/nam_micro_lookup.dat ./eta_micro_lookup.dat

export PARMnam=$homedir/parm

# copy flat files instead
cp ${svndir}/parm/postxconfig-NT-NMM.txt ./postxconfig-NT.txt
cp ${svndir}/parm/params_grib2_tbl_new params_grib2_tbl_new

${APRUN} ${POSTGPEXEC} < itag > outpost_nems_${NEWDATE}

mv BGDAWP${fhr}.tm00 BGDAWP${fhr}.tm00.Grib2
mv BGRD3D${fhr}.tm00 BGRD3D${fhr}.tm00.Grib2
mv BGRDSF${fhr}.tm00 BGRDSF${fhr}.tm00.Grib2

# operational NMMB post processing generates 3 files
filelist="BGDAWP${fhr}.tm00.Grib2 \
          BGRD3D${fhr}.tm00.Grib2 \
          BGRDSF${fhr}.tm00.Grib2"

for file in $filelist; do
export filein2=$file
ls -l ${filein2}
export err=$?

if [ $err = "0" ] ; then

 # operational NMMB post processing generates 3 files, start with BGDAWP first
 # use cmp to see if new pgb files are identical to the control one
 cmp ${filein2} $homedir/data_out/nmmb/${filein2}.${machine}

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
  $cmp_grib2_grib2 $homedir/data_out/nmmb/${filein2}.${machine} ${filein2} > ${filein2}.diff
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
