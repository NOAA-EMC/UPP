#!/bin/bash

#SBATCH -o out.post.aigfs
#SBATCH -e out.post.aigfs
#SBATCH -J aigfs_test 
#SBATCH -t @[WTIME]
#SBATCH -q @[QUEUE]
#SBATCH -A @[accnr]
#SBATCH @[EXCLUSIVE]
#SBATCH @[N_TASKS]
#SBATCH @[TASKS_PER_NODE]
#SBATCH @[NODES] @[N_TASKS_PER_NODE]

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
module use ${svndir}/modulefiles
module load $(echo "${machine}" | tr '[:upper:]' '[:lower:]')_${compiler}
module load wgrib2/3.6.0
module load prod_util/2.1.1
module list

@[STACK_SIZE]

msg="Starting aigfs test"
postmsg "$logfile" "$msg"

export POSTGPEXEC=${svndir}/exec/upp.x

# specify forecast start time and hour for running your post job
export startdate=2026050518
export fhr=084
export cyc=`echo $startdate | cut -c9-10`

# specify your running and output directory
export DATA=$rundir/aigfs_${startdate}
rm -rf $DATA; mkdir -p $DATA
cd $DATA

#export NEWDATE=`${NDATE} +${fhr} $startdate` 
NEWDATE=$startdate
export YY=`echo $NEWDATE | cut -c1-4`
export MM=`echo $NEWDATE | cut -c5-6`
export DD=`echo $NEWDATE | cut -c7-8`
export HH=`echo $NEWDATE | cut -c9-10`

cat > itag <<EOF
&model_inputs
fileName='$homedir/data_in/aigfs/aigfs.t${cyc}z.f${fhr}.nc'
IOFORM='netcdf'
grib='grib2'
DateStr='${YY}-${MM}-${DD}_${HH}:00:00'
MODELNAME='AIGFS'
fhr=${fhr}
/
&NAMPGB
KPO=13,PO=50.,100.,150.,200.,250.,300.,400.,500.,600.,700.,850.,925.,1000.,
/
EOF

# copy fix data
cp ${svndir}/fix/nam_micro_lookup.dat ./eta_micro_lookup.dat
cp ${svndir}/parm/aigfs/postxconfig-NT-aigfs.txt ./postxconfig-NT.txt
cp ${svndir}/parm/params_grib2_tbl_new ./params_grib2_tbl_new

# Run the UPP
${APRUN} ${POSTGPEXEC} < itag > outpost_aigfs_${NEWDATE}

################################################
# Compare with baseline data
################################################
fhr=$((10#$fhr))
FH3=$(printf "%03d" "$fhr")
FH2=$(printf "%02d" "$fhr")
mv PRS.GrbF${FH2} aigfs.t${cyc}z.pres.f${FH3}.grib2
mv SFC.GrbF${FH2} aigfs.t${cyc}z.sfc.f${FH3}.grib2

# AIGFS post processing generates 2 file
filelist="aigfs.t${cyc}z.pres.f${FH3}.grib2 \
          aigfs.t${cyc}z.sfc.f${FH3}.grib2"

for file in $filelist; do
export filein2=$file
ls -l ${filein2}
export err=$?

if [ $err = "0" ] ; then

 # use cmp to see if new pgb files are identical to the control one
 cmp ${filein2} $homedir/data_out_$compiler/aigfs/${filein2}.${machine}

 # if not bit-identical, use cmp_grib2_grib2 to compare each grib record
 export err1=$?
 if [ $err1 -eq 0 ] ; then
  msg="aigfs test: your new post executable generates bit-identical ${filein2} as the develop branch"
  echo $msg
 else
  msg="aigfs test: your new post executable did not generate bit-identical ${filein2} as the develop branch"
  echo $msg
  echo " start comparing each grib record and write the comparison result to *diff files"
  echo " check these *diff files to make sure your new post only change variables which you intend to change"
  $cmp_grib2_grib2 $homedir/data_out_$compiler/aigfs/${filein2}.${machine} ${filein2} > ${filein2}.diff
 fi
else
 msg="aigfs test: post failed using your new post executable to generate ${filein2}"
 echo $msg 2>&1 | tee -a TEST_ERROR
fi

postmsg "$logfile" "$msg"
done

echo "PROGRAM IS COMPLETE!!!!!" 2>&1 | tee SUCCESS
msg="Ending aigfs test"
postmsg "$logfile" "$msg"
