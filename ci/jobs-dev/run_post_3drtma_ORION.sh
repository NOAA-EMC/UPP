#!/bin/sh 
 
#SBATCH -o out.post.rtma
#SBATCH -e out.post.rtma
#SBATCH -J rtma_test
#SBATCH -t 00:30:00
#SBATCH -N 8 --ntasks-per-node=12
#SBATCH -q batch
#SBATCH -A nems

set -x

# specify computation resource
export threads=1
export OMP_NUM_THREADS=$threads
export MP_LABELIO=yes
export APRUN="srun"

echo "starting time"
date

######################################################################
# Purpose: to run RAP post processing
######################################################################

# EXPORT list here

module use /apps/contrib/spack-stack/spack-stack-1.8.0/envs/ue-intel-2021.9.0/install/modulefiles/Core
module load stack-intel/2021.9.0
module load stack-intel-oneapi-mpi/2021.9.0
module load libpng/1.6.37
module load jasper/2.0.32
module load prod_util/2.1.1
module load crtm/2.4.0.1
module list

ulimit -s unlimited
export WGRIB2=wgrib2
export COMROOT=$rundir

msg="Starting rtma test"
postmsg "$logfile" "$msg"

export cmp_grib2_grib2=/home/wmeng/bin/cmp_grib2_grib2_new
export POSTGPEXEC=${svndir}/exec/upp.x

# CALL executable job script here

# specify your running and output directory
export startdate=2023040400
export fhr=000
export tmmark=tm00
export DATA=$rundir/rtma_${startdate}

export NEWDATE=$startdate

export YY=`echo ${NEWDATE} | cut -c1-4`
export MM=`echo ${NEWDATE} | cut -c5-6`
export DD=`echo ${NEWDATE} | cut -c7-8`
export HH=`echo ${NEWDATE} | cut -c9-10`
export min=00

rm -rf $DATA; mkdir -p $DATA
cd $DATA

cat > itag <<EOF
&model_inputs
fileName='$homedir/data_in/3drtma/dynf${fhr}.nc'
IOFORM='netcdf'
grib='grib2'
DateStr='${YY}-${MM}-${DD}_${HH}:${min}:00'
MODELNAME='FV3R'
SUBMODELNAME='RTMA'
fileNameFlux='$homedir/data_in/3drtma/phyf${fhr}.nc'
/
&NAMPGB
KPO=47,PO=1000.,975.,950.,925.,900.,875.,850.,825.,800.,775.,750.,725.,700.,675.,650.,625.,600.,575.,550.,525.,500.,475.,450.,425.,400.,375.,350.,325.,300.,275.,250.,225.,200.,175.,150.,125.,100.,70.,50.,30.,20.,10.,7.,5.,3.,2.,1.,
/
EOF

#copy fix data
cp ${svndir}/parm/params_grib2_tbl_new params_grib2_tbl_new
cp ${svndir}/parm/postxconfig-NT-rrfs.txt postxconfig-NT.txt
cp ${svndir}/fix/nam_micro_lookup.dat eta_micro_lookup.dat

#get crtm fix file
for what in "amsre_aqua" "imgr_g11" "imgr_g12" "imgr_g13" \
    "imgr_g15" "imgr_mt1r" "imgr_mt2" "seviri_m10" \
    "ssmi_f13" "ssmi_f14" "ssmi_f15" "ssmis_f16" \
    "ssmis_f17" "ssmis_f18" "ssmis_f19" "ssmis_f20" \
    "tmi_trmm" "v.seviri_m10" "imgr_insat3d" "abi_gr" \
    "ahi_himawari8" ; do
    ln -s "${CRTM_FIX}/${what}.TauCoeff.bin" .
    ln -s "${CRTM_FIX}/${what}.SpcCoeff.bin" .
done

for what in 'Aerosol' 'Cloud' ; do
    ln -s "${CRTM_FIX}/${what}Coeff.bin" .
done

for what in  ${CRTM_FIX}/*Emis* ; do
   ln -s $what .
done

#copy xml
cp ${svndir}/parm/params_grib2_tbl_new params_grib2_tbl_new
cp ${svndir}/parm/postxconfig-NT-rrfs.txt postxconfig-NT.txt
cp ${svndir}/fix/nam_micro_lookup.dat eta_micro_lookup.dat

${APRUN} ${POSTGPEXEC} < itag > wrfpost2.out

# operational rtma post processing generates 3 files
filelist="NATLEV00.tm00 \
          PRSLEV00.tm00 \
          IFIFIP00.tm00"

for file in $filelist; do
export filein2=$file
ls -l ${filein2}
export err=$?

if [ $err = "0" ] ; then

 # operational rtma post processing generates 3 files, start with BGDAWP first
 # use cmp to see if new pgb files are identical to the control one
 cmp ${filein2} $homedir/data_out/3drtma/${filein2}.${machine}

 # if not bit-identical, use cmp_grib2_grib2 to compare each grib record
 export err1=$?
 if [ $err1 -eq 0 ] ; then
  msg="rtma test: your new post executable generates bit-identical ${filein2} as the trunk"
  echo $msg
 else
  msg="rtma test: your new post executable did not generate bit-identical ${filein2} as the trunk"
  echo $msg
  $cmp_grib2_grib2 $homedir/data_out/3drtma/${filein2}.${machine} ${filein2} > ${filein2}.diff
 fi


else

 msg="rtma test: post failed using your new post executable to generate ${filein2}"
 echo $msg 2>&1 | tee -a TEST_ERROR

fi
postmsg "$logfile" "$msg"
done

echo "PROGRAM IS COMPLETE!!!!!" 2>&1 | tee SUCCESS
msg="Ending rtma test"
postmsg "$logfile" "$msg"


