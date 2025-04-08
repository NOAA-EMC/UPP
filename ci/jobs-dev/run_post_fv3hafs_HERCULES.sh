#!/bin/sh

#SBATCH -o out.fv3hafs
#SBATCH -e out.fv3hafs
#SBATCH -J fv3hafs_test 
#SBATCH -t 00:20:00
#SBATCH -N 5 --ntasks-per-node=12
#SBATCH -q batch
#SBATCH -A nems
#SBATCH --exclusive

set -x

# specify computation resource
export threads=1
#export MP_LABELIO=yes
export OMP_NUM_THREADS=$threads
export APRUN="srun"

############################################
# Loading module
############################################
module use /apps/contrib/spack-stack/spack-stack-1.8.0/envs/ue-intel-2021.9.0/install/modulefiles/Core
module load stack-intel/2021.9.0
module load stack-intel-oneapi-mpi/2021.9.0
module load libpng/1.6.37
module load jasper/2.0.32
module load prod_util/2.1.1
module load crtm/2.4.0.1
module list

#export WGRIB2=wgrib2
#export COMROOT=$rundir
#export CRTM_FIX=/apps/contrib/NCEPLIBS/orion/fix/crtm_v2.3.0

#ulimit -s1900000000
ulimit -s unlimited

msg="Starting fv3hafs test"
postmsg "$logfile" "$msg"

export cmp_grib2_grib2=/home/wmeng/bin/cmp_grib2_grib2_new
export POSTGPEXEC=${svndir}/exec/upp.x     

# specify forecast start time and hour for running your post job
export startdate=2022092800
export fhr=009
export CC=`echo $startdate | cut -c9-10`

# specify your running and output directory
export DATA=$rundir/fv3hafs_${startdate}
export tmmark=tm00
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

rm -f fort.*

cp ${svndir}/fix/nam_micro_lookup.dat ./eta_micro_lookup.dat

# copy flat files instead
cp ${svndir}/parm/postxconfig-NT-hafs_nosat.txt ./postxconfig-NT.txt

cp ${svndir}/parm/params_grib2_tbl_new ./params_grib2_tbl_new

${APRUN} ${POSTGPEXEC} < itag > outpost_nems_${NEWDATE}

fhr=`expr $fhr + 0`
fhr2=`printf "%02d" $fhr`

filelist="HURPRS${fhr2}.tm00"

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
  msg="fv3hafs test: your new post executable generates bit-identical ${filein2} as the trunk"
  echo $msg
 else
  msg="fv3hafs test: your new post executable did not generate bit-identical ${filein2} as the trunk"
  echo $msg
  echo " start comparing each grib record and write the comparison result to *diff files"
  echo " check these *diff files to make sure your new post only change variables which you intend to change"
  $cmp_grib2_grib2 $homedir/data_out/hafs/${filein2}.${machine} ${filein2} > ${filein2}.diff
 fi

else

 msg="fv3hafs test: post failed using your new post executable to generate ${filein2}"
 echo $msg 2>&1 | tee -a TEST_ERROR

fi
postmsg "$logfile" "$msg"
done

echo "PROGRAM IS COMPLETE!!!!!" 2>&1 | tee SUCCESS
msg="Ending fv3hafs test"
postmsg "$logfile" "$msg"
