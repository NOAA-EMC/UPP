#!/bin/bash

#PBS -o out.post.dafs
#PBS -e out.post.dafs
#PBS -N dafs_test
#PBS -l walltime=00:30:00
#PBS -q debug
#PBS -A GFS-DEV
#PBS -l place=vscatter,select=2:ncpus=60:mem=300GB
#PBS -V

set -x

# specify computation resources
export threads=1
export OMP_NUM_THREADS=$threads
export APRUN="mpiexec -l -n 120 -ppn 60"

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

msg="Starting dafs test"
postmsg "$logfile" "$msg"

export POSTGPEXEC=${svndir}/exec/upp.x

# specify forecast start time and hour for running your post job
export startdate=2025061712
export fhr=08

# specify your running and output directory
export DATA=$rundir/dafs_${startdate}
rm -rf $DATA; mkdir -p $DATA
cd $DATA

export NEWDATE=`${NDATE} +${fhr} $startdate`
export YY=`echo $NEWDATE | cut -c1-4`
export MM=`echo $NEWDATE | cut -c5-6`
export DD=`echo $NEWDATE | cut -c7-8`
export HH=`echo $NEWDATE | cut -c9-10`

cat > itag <<EOF
&model_inputs
fileName='$homedir/data_in/dafs/wrfout_d01_${YY}-${MM}-${DD}_${HH}_00_00'
IOFORM='netcdf'
grib='grib2'
DateStr='${YY}-${MM}-${DD}_${HH}:00:00'
MODELNAME='RAPR'
SUBMODELNAME='RAPR'
/
&NAMPGB
KPO=47,PO=2.,5.,7.,10.,20.,30.,50.,70.,75.,100.,125.,150.,175.,200.,225.,250.,275.,300.,325.,350.,375.,400.,425.,450.,475.,500.,525.,550.,575.,600.,625.,650.,675.,700.,725.,750.,775.,800.,825.,850.,875.,900.,925.,950.,975.,1000.,1013.2,gtg_on=.true.,
/
EOF

# copy fix data
cp ${svndir}/fix/nam_micro_lookup.dat                      ./eta_micro_lookup.dat
cp ${svndir}/parm/dafs/postxconfig-NT-hrrr_dafs.txt        ./postxconfig-NT.txt
cp ${svndir}/parm/params_grib2_tbl_new                     ./params_grib2_tbl_new
cp ${svndir}/sorc/ncep_post.fd/post_gtg.fd/gtg.config.hrrr ./gtg.config.hrrr
cp ${svndir}/sorc/ncep_post.fd/post_gtg.fd/gtg.input.hrrr  ./gtg.input.hrrr

# Run the UPP
${APRUN} ${POSTGPEXEC} < itag > outpost_dafs_${NEWDATE}

################################################
# Compare with baseline data
################################################
fhr=`expr $fhr + 0`
fhr2=`printf "%02d" $fhr`

# DAFS post processing generates 2 files
filelist="AVIATION.GrbF${fhr2} \
          IFIFIP.GrbF${fhr2}"

for file in $filelist; do
export filein2=$file
ls -l ${filein2}
export err=$?

if [ $err = "0" ] ; then

 # use cmp to see if new pgb files are identical to the control one
 cmp ${filein2} $homedir/data_out/dafs/${filein2}.${machine}

 # if not bit-identical, use cmp_grib2_grib2 to compare each grib record
 export err1=$?
 if [ $err1 -eq 0 ] ; then
  msg="dafs test: your new post executable generates bit-identical ${filein2} as the develop branch"
  echo $msg
 else
  msg="dafs test: your new post executable did not generate bit-identical ${filein2} as the develop branch"
  echo $msg
  echo " start comparing each grib record and write the comparison result to *diff files"
  echo " check these *diff files to make sure your new post only change variables which you intend to change"
  $cmp_grib2_grib2 $homedir/data_out/dafs/${filein2}.${machine} ${filein2} > ${filein2}.diff
 fi
else
 msg="dafs test: post failed using your new post executable to generate ${filein2}"
 echo $msg 2>&1 | tee -a TEST_ERROR
fi

postmsg "$logfile" "$msg"
done

echo "PROGRAM IS COMPLETE!!!!!" 2>&1 | tee SUCCESS
msg="Ending rrfs_ifi test"
postmsg "$logfile" "$msg"
