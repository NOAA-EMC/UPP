#!/bin/bash

#PBS -o out.post.aqm
#PBS -e out.post.aqm
#PBS -N aqm_test
#PBS -l walltime=00:30:00
#PBS -q debug
#PBS -A GFS-DEV
#PBS -l place=vscatter,select=2:ncpus=48
#PBS -V

set -x

# specify computation resources
export threads=1
export OMP_NUM_THREADS=$threads
export APRUN="mpiexec -l -n 96 -ppn 48"

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

msg="Starting aqm test"
postmsg "$logfile" "$msg"

export POSTGPEXEC=${svndir}/exec/upp.x

# specify forecast start time and hour for running your post job
export startdate=2025103106
export fhr=008
export cyc=`echo $startdate | cut -c9-10`
export tmmark=tm00

# specify your running and output directory
export DATA=$rundir/aqm_${startdate}
rm -rf $DATA; mkdir -p $DATA
cd $DATA

export NEWDATE=`${NDATE} +${fhr} $startdate` 
export YY=`echo $NEWDATE | cut -c1-4`
export MM=`echo $NEWDATE | cut -c5-6`
export DD=`echo $NEWDATE | cut -c7-8`
export HH=`echo $NEWDATE | cut -c9-10`

cat > itag <<EOF
&model_inputs
fileName='$homedir/data_in/aqm/aqm.t${cyc}z.dyn.f${fhr}.nc'
IOFORM='netcdf'
grib='grib2'
DateStr='${YY}-${MM}-${DD}_${HH}:00:00'
MODELNAME='FV3R'
fileNameFlux='$homedir/data_in/aqm/aqm.t${cyc}z.phy.f${fhr}.nc'
/
&NAMPGB
KPO=47,PO=1000.,975.,950.,925.,900.,875.,850.,825.,800.,775.,750.,725.,700.,675.,650.,625.,600.,575.,550.,525.,500.,475.,450.,425.,400.,375.,350.,325.,300.,275.,250.,225.,200.,175.,150.,125.,100.,70.,50.,30.,20.,10.,7.,5.,3.,2.,1.,aqf_on=.true.,
/
EOF

# copy fix data
cp ${svndir}/fix/nam_micro_lookup.dat ./eta_micro_lookup.dat
cp ${svndir}/parm/aqm/postxconfig-NT-AQM.txt ./postxconfig-NT.txt
cp ${svndir}/parm/params_grib2_tbl_new ./params_grib2_tbl_new

# Run the UPP
${APRUN} ${POSTGPEXEC} < itag > outpost_aqm_${NEWDATE}

################################################
# Compare with baseline data
################################################
fhr=`expr $fhr + 0`
fhr2=`printf "%02d" $fhr`

# AQM post processing generates 1 file
filelist="CMAQ${fhr2}.${tmmark}"

for file in $filelist; do
export filein2=$file
ls -l ${filein2}
export err=$?

if [ $err = "0" ] ; then
 # use cmp to see if new pgb files are identical to the control one
 cmp ${filein2} $homedir/data_out/aqm/${filein2}.${machine}

 # if not bit-identical, use cmp_grib2_grib2 to compare each grib record
 export err1=$?
 if [ $err1 -eq 0 ] ; then
  msg="aqm test: your new post executable generates bit-identical ${filein2} as the develop branch"
  echo $msg
 else
  msg="aqm test: your new post executable did not generate bit-identical ${filein2} as the develop branch"
  echo $msg
  echo " start comparing each grib record and write the comparison result to *diff files"
  echo " check these *diff files to make sure your new post only change variables which you intend to change"
  $cmp_grib2_grib2 $homedir/data_out/aqm/${filein2}.${machine} ${filein2} > ${filein2}.diff
 fi
else
 msg="aqm test: post failed using your new post executable to generate ${filein2}"
 echo $msg 2>&1 | tee -a TEST_ERROR
fi

postmsg "$logfile" "$msg"
done

echo "PROGRAM IS COMPLETE!!!!!" 2>&1 | tee SUCCESS
msg="Ending aqm test"
postmsg "$logfile" "$msg"
