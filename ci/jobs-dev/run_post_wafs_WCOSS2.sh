#!/bin/bash

#PBS -o out.post.wafs
#PBS -e out.post.wafs
#PBS -N wafs_test
#PBS -l walltime=00:30:00
#PBS -q debug
#PBS -A GFS-DEV
#PBS -l place=vscatter,select=2:ncpus=48:mem=300GB
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

msg="Starting wafs test"
postmsg "$logfile" "$msg"

export POSTGPEXEC=${svndir}/exec/upp.x

# specify forecast start time and hour for running your post job
export startdate=2025012600
export fhr=006
export cyc=`echo $startdate |cut -c9-10`

# specify your running and output directory
export DATA=$rundir/wafs_${startdate}
rm -rf $DATA; mkdir -p $DATA
cd $DATA

export NEWDATE=`${NDATE} +${fhr} $startdate`
export YY=`echo $NEWDATE | cut -c1-4`
export MM=`echo $NEWDATE | cut -c5-6`
export DD=`echo $NEWDATE | cut -c7-8`
export HH=`echo $NEWDATE | cut -c9-10`

cat > itag <<EOF
&model_inputs
fileName='$homedir/data_in/wafs/gfs.t${cyc}z.atmf${fhr}.nc'
IOFORM='netcdf'
grib='grib2'
DateStr='${YY}-${MM}-${DD}_${HH}:00:00'
MODELNAME='GFS'
SUBMODELNAME='GFS'
fileNameFlux='$homedir/data_in/wafs/gfs.t${cyc}z.sfcf${fhr}.nc'
/
&NAMPGB
KPO=60,PO=97720.,94210.,90810.,87510.,84310.,81200.,78190.,75260.,72430.,69680.,67020.,64440.,61940.,59520.,57180.,54920.,52720.,50600.,48550.,46560.,44650.,42790.,41000.,39270.,37600.,35990.,34430.,32930.,31490.,30090.,28740.,27450.,26200.,25000.,23840.,22730.,21660.,20650.,19680.,18750.,17870.,17040.,16240.,15470.,14750.,14060.,13400.,12770.,12170.,11600.,11050.,10530.,10040.,9570.,9120.,8700.,8280.,7900.,7520.,7170.,gtg_on=.true.,popascal=.true.,
/
EOF

# copy fix data
cp ${svndir}/fix/nam_micro_lookup.dat                     ./eta_micro_lookup.dat
cp ${svndir}/parm/wafs/postxconfig-NT-gfs-wafs.txt        ./postxconfig-NT.txt
cp ${svndir}/parm/params_grib2_tbl_new                    ./params_grib2_tbl_new
cp ${svndir}/sorc/ncep_post.fd/post_gtg.fd/gtg.config.gfs .
cp ${svndir}/sorc/ncep_post.fd/post_gtg.fd/gtg.input.gfs  .
cp ${svndir}/sorc/ncep_post.fd/post_gtg.fd/imprintings.gtg_gfs.txt .

# Run the UPP
${APRUN} ${POSTGPEXEC} < itag > outpost_wafs_${NEWDATE}

################################################
# Compare with baseline data
################################################
fhr=`expr $fhr + 0`
fhr2=`printf "%02d" $fhr`

# WAFS post processing generates 1 file
filelist="GFSPRS.GrbF${fhr2}"

for file in $filelist; do
export filein2=$file
ls -l ${filein2}
export err=$?

if [ $err = "0" ] ; then

 # use cmp to see if new pgb files are identical to the control one
 cmp ${filein2} $homedir/data_out/wafs/${filein2}.${machine}

 # if not bit-identical, use cmp_grib2_grib2 to compare each grib record
 export err1=$?
 if [ $err1 -eq 0 ] ; then
  msg="wafs test: your new post executable generates bit-identical ${filein2} as the develop branch"
  echo $msg
 else
  msg="wafs test: your new post executable did not generate bit-identical ${filein2} as the develop branch"
  echo $msg
  echo " start comparing each grib record and write the comparison result to *diff files"
  echo " check these *diff files to make sure your new post only change variables which you intend to change"
  $cmp_grib2_grib2 $homedir/data_out/wafs/${filein2}.${machine} ${filein2} > ${filein2}.diff
 fi
else
 msg="wafs test: post failed using your new post executable to generate ${filein2}"
 echo $msg 2>&1 | tee -a TEST_ERROR
fi

postmsg "$logfile" "$msg"
done

echo "PROGRAM IS COMPLETE!!!!!" 2>&1 | tee SUCCESS
msg="Ending wafs test"
postmsg "$logfile" "$msg"
