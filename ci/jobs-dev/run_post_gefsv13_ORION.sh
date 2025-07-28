#!/bin/bash

#SBATCH -o out.post.gefsv13
#SBATCH -e out.post.gefsv13
#SBATCH -J gefsv13_test 
#SBATCH -t 00:30:00
#SBATCH -N 3 --ntasks-per-node=12
#SBATCH -q batch
#SBATCH -A nems

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
module load orion_$compiler
module load wgrib2/3.6.0
module load prod_util/2.1.1
module list

ulimit -s unlimited

msg="Starting gefsv13 test"
postmsg "$logfile" "$msg"

export POSTGPEXEC=${svndir}/exec/upp.x

# specify forecast start time and hour for running your post job
export startdate=2023080300
export fhr=009
export cyc=`echo $startdate | cut -c9-10`

# specify your running and output directory
export DATA=$rundir/gefsv13_${startdate}
rm -rf $DATA; mkdir -p $DATA
cd $DATA

export NEWDATE=`${NDATE} +${fhr} $startdate` 
export YY=`echo $NEWDATE | cut -c1-4`
export MM=`echo $NEWDATE | cut -c5-6`
export DD=`echo $NEWDATE | cut -c7-8`
export HH=`echo $NEWDATE | cut -c9-10`

cat > itag <<EOF
&model_inputs
fileName='$homedir/data_in/gefsv13/gefs.t${cyc}z.atmf${fhr}.nc'
IOFORM='netcdf'
grib='grib2'
DateStr='${YY}-${MM}-${DD}_${HH}:00:00'
MODELNAME='FV3R'
fileNameFlux='$homedir/data_in/gefsv13/gefs.t${cyc}z.sfcf${fhr}.nc'
/
&NAMPGB
KPO=50,PO=1000.,975.,950.,925.,900.,875.,850.,825.,800.,775.,750.,725.,700.,675.,650.,625.,600.,575.,550.,525.,500.,475.,450.,425.,400.,375.,350.,325.,300.,275.,250.,225.,200.,175.,150.,125.,100.,70.,50.,40.,30.,20.,15.,10.,7.,5.,3.,2.,1.,0.4,
/
EOF

export e1=3
export e2=01
export e3=30

# copy fix data
cp ${svndir}/fix/nam_micro_lookup.dat ./eta_micro_lookup.dat
cp ${svndir}/parm/gefs/postxconfig-NT-gefs.txt ./postxconfig-NT.txt
cp ${svndir}/parm/params_grib2_tbl_new ./params_grib2_tbl_new

# Run the UPP
${APRUN} ${POSTGPEXEC} < itag > outpost_gefsv13_${NEWDATE}

################################################
# Compare with baseline data
################################################
fhr=$((10#$fhr))
FH3=$(printf "%03d" "$fhr")
FH2=$(printf "%02d" "$fhr")
mv GFSPRS.GrbF${FH2} gefs.t${cyc}z.master.grb2f${FH3}

# GEFSv13 post processing generates 1 file
filelist="gefs.t${cyc}z.master.grb2f${FH3}"

for file in $filelist; do
export filein2=$file
ls -l ${filein2}
export err=$?

if [ $err = "0" ] ; then
 # use cmp to see if new pgb files are identical to the control one
 cmp ${filein2} $homedir/data_out_$compiler/gefsv13/${filein2}.${machine}

 # if not bit-identical, use cmp_grib2_grib2 to compare each grib record
 export err1=$?
 if [ $err1 -eq 0 ] ; then
  msg="gefsv13 test: your new post executable generates bit-identical ${filein2} as the develop branch"
  echo $msg
 else
  msg="gefsv13 test: your new post executable did not generate bit-identical ${filein2} as the develop branch"
  echo $msg
  echo " start comparing each grib record and write the comparison result to *diff files"
  echo " check these *diff files to make sure your new post only change variables which you intend to change"
  $cmp_grib2_grib2 $homedir/data_out_$compiler/gefsv13/${filein2}.${machine} ${filein2} > ${filein2}.diff
 fi
else
 msg="gefsv13 test: post failed using your new post executable to generate ${filein2}"
 echo $msg 2>&1 | tee -a TEST_ERROR
fi

postmsg "$logfile" "$msg"
done

echo "PROGRAM IS COMPLETE!!!!!" 2>&1 | tee SUCCESS
msg="Ending gefsv13 test"
postmsg "$logfile" "$msg"
