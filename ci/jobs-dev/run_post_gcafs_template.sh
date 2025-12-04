#!/bin/bash

#SBATCH -o out.post.gcafs
#SBATCH -e out.post.gcafs
#SBATCH -J gcafs_test 
#SBATCH -t @[WTIME]
#SBATCH -q @[QUEUE]
#SBATCH -A @[accnr]
#SBATCH @[EXCLUSIVE]
#SBATCH -N @[NODES] --ntasks-per-node=@[N_TASKS_PER_NODE]

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

msg="Starting gcafs test"
postmsg "$logfile" "$msg"

export POSTGPEXEC=${svndir}/exec/upp.x

# specify forecast start time and hour for running your post job
export startdate=2025110500
export fhr=060
export cyc=`echo $startdate | cut -c9-10`

# specify your running and output directory
export DATA=$rundir/gcafs_${startdate}
rm -rf $DATA; mkdir -p $DATA
cd $DATA

export NEWDATE=`${NDATE} +${fhr} $startdate` 
export YY=`echo $NEWDATE | cut -c1-4`
export MM=`echo $NEWDATE | cut -c5-6`
export DD=`echo $NEWDATE | cut -c7-8`
export HH=`echo $NEWDATE | cut -c9-10`

cat > itag <<EOF
&model_inputs
fileName='$homedir/data_in/gcafs/gcafs.t${cyc}z.atmf${fhr}.nc'
IOFORM='netcdf'
grib='grib2'
DateStr='${YY}-${MM}-${DD}_${HH}:00:00'
MODELNAME='FV3R'
fileNameFlux='$homedir/data_in/gcafs/gcafs.t${cyc}z.sfcf${fhr}.nc'
/
&NAMPGB
KPO=57,PO=1000.,975.,950.,925.,900.,875.,850.,825.,800.,775.,750.,725.,700.,675.,650.,625.,600.,575.,550.,525.,500.,475.,450.,425.,400.,375.,350.,325.,300.,275.,250.,225.,200.,175.,150.,125.,100.,70.,50.,40.,30.,20.,15.,10.,7.,5.,3.,2.,1.,0.7,0.4,0.2,0.1,0.07,0.04,0.02,0.01,nasa_on=.true.,
/
EOF


# copy fix data
cp ${svndir}/fix/nam_micro_lookup.dat ./eta_micro_lookup.dat
cp ${svndir}/parm/gcafs/postxconfig-NT-gcafs.txt ./postxconfig-NT.txt
cp ${svndir}/parm/params_grib2_tbl_new ./params_grib2_tbl_new

cp ${svndir}/fix/chem/optics_luts_DUST_nasa.dat .
cp ${svndir}/fix/chem/optics_luts_SALT_nasa.dat .
cp ${svndir}/fix/chem/optics_luts_SOOT_nasa.dat .
cp ${svndir}/fix/chem/optics_luts_SUSO_nasa.dat .
cp ${svndir}/fix/chem/optics_luts_WASO_nasa.dat .
cp ${svndir}/fix/chem/optics_luts_NITR_nasa.dat .

# Run the UPP
${APRUN} ${POSTGPEXEC} < itag > outpost_gcafs_${NEWDATE}

################################################
# Compare with baseline data
################################################
fhr=$((10#$fhr))
FH3=$(printf "%03d" "$fhr")
FH2=$(printf "%02d" "$fhr")

# GCAFS post processing generates 2 file
filelist="GFSPRS.GrbF${FH2} \
          GFSFLX.GrbF${FH2}"

for file in $filelist; do
export filein2=$file
ls -l ${filein2}
export err=$?

if [ $err = "0" ] ; then
 # use cmp to see if new pgb files are identical to the control one
 cmp ${filein2} $homedir/data_out_$compiler/gcafs/${filein2}.${machine}

 # if not bit-identical, use cmp_grib2_grib2 to compare each grib record
 export err1=$?
 if [ $err1 -eq 0 ] ; then
  msg="gcafs test: your new post executable generates bit-identical ${filein2} as the develop branch"
  echo $msg
 else
  msg="gcafs test: your new post executable did not generate bit-identical ${filein2} as the develop branch"
  echo $msg
  echo " start comparing each grib record and write the comparison result to *diff files"
  echo " check these *diff files to make sure your new post only change variables which you intend to change"
  $cmp_grib2_grib2 $homedir/data_out_$compiler/gcafs/${filein2}.${machine} ${filein2} > ${filein2}.diff
 fi
else
 msg="gcafs test: post failed using your new post executable to generate ${filein2}"
 echo $msg 2>&1 | tee -a TEST_ERROR
fi

postmsg "$logfile" "$msg"
done

echo "PROGRAM IS COMPLETE!!!!!" 2>&1 | tee SUCCESS
msg="Ending gcafs test"
postmsg "$logfile" "$msg"
