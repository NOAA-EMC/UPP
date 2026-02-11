#!/bin/bash

#PBS -o out.post.hafs
#PBS -e out.post.hafs
#PBS -N hafs_test
#PBS -l walltime=00:30:00
#PBS -q debug
#PBS -A GFS-DEV
#PBS -l place=vscatter,select=3:ncpus=24
#PBS -V

set -x

# specify computation resources
export threads=1
export OMP_NUM_THREADS=$threads
export APRUN="mpiexec -l -n 72 -ppn 24"

echo "starting time"
date

############################################
# Loading modules
############################################
module reset
module use ${svndir}/modulefiles
module load wcoss2_intel
module load cray-pals/1.0.12
module load libjpeg/9c
module load prod_util/2.0.14
module load wgrib2/2.0.8
module list

msg="Starting hafs test"
postmsg "$logfile" "$msg"

export POSTGPEXEC=${svndir}/exec/upp.x

# specify forecast start time and hour for running your post job
export startdate=2022092800
export fhr=009
export tmmark=tm00

# specify your running and output directory
export DATA=$rundir/hafs_${startdate}
rm -rf $DATA; mkdir -p $DATA
cd $DATA

export NEWDATE=`${NDATE} +${fhr} $startdate` 
export YY=`echo $NEWDATE | cut -c1-4`
export MM=`echo $NEWDATE | cut -c5-6`
export DD=`echo $NEWDATE | cut -c7-8`
export HH=`echo $NEWDATE | cut -c9-10`

cat > itag <<EOF
&model_inputs
fileName='$homedir/data_in/hafs/atmf${fhr}.nc'
IOFORM='netcdf'
grib='grib2'
DateStr='${YY}-${MM}-${DD}_${HH}:00:00'
MODELNAME='FV3R'
fileNameFlux='$homedir/data_in/hafs/sfcf${fhr}.nc'
/
&NAMPGB
KPO=47,PO=1000.,975.,950.,925.,900.,875.,850.,825.,800.,775.,750.,725.,700.,675.,650.,625.,600.,575.,550.,525.,500.,475.,450.,425.,400.,375.,350.,325.,300.,275.,250.,225.,200.,175.,150.,125.,100.,70.,50.,30.,20.,10.,7.,5.,3.,2.,1.,
/
EOF

# copy fix data
cp ${svndir}/fix/nam_micro_lookup.dat ./eta_micro_lookup.dat
cp ${svndir}/parm/hafs/postxconfig-NT-hafs_nosat.txt ./postxconfig-NT.txt
cp ${svndir}/parm/params_grib2_tbl_new ./params_grib2_tbl_new

# Run the UPP
${APRUN} ${POSTGPEXEC} < itag > outpost_hafs_${NEWDATE}

################################################
# Compare with baseline data
################################################
fhr=`expr $fhr + 0`
fhr2=`printf "%02d" $fhr`

# HAFS post processing generates 1 file
filelist="HURPRS${fhr2}.${tmmark}"

for file in $filelist; do
export filein2=$file
ls -l ${filein2}
export err=$?

if [ $err = "0" ] ; then
 # use cmp to see if new pgb files are identical to the control one
 cmp ${filein2} $homedir/data_out/hafs/${filein2}.${machine}

 # if not bit-identical, use cmp_grib2_grib2 to compare each grib record
 export err1=$?
 if [ $err1 -eq 0 ] ; then
  msg="hafs test: your new post executable generates bit-identical ${filein2} as the develop branch"
  echo $msg
 else
  msg="hafs test: your new post executable did not generate bit-identical ${filein2} as the develop branch"
  echo $msg
  echo " start comparing each grib record and write the comparison result to *diff files"
  echo " check these *diff files to make sure your new post only change variables which you intend to change"
  $cmp_grib2_grib2 $homedir/data_out/hafs/${filein2}.${machine} ${filein2} > ${filein2}.diff
 fi
else
 msg="hafs test: post failed using your new post executable to generate ${filein2}"
 echo $msg 2>&1 | tee -a TEST_ERROR
fi

postmsg "$logfile" "$msg"
done

echo "PROGRAM IS COMPLETE!!!!!" 2>&1 | tee SUCCESS
msg="Ending hafs test"
postmsg "$logfile" "$msg"
