#!/bin/bash 
 
#PBS -o out.post.3drtma
#PBS -e out.post.3drtma
#PBS -N 3drtma.test
#PBS -l walltime=00:30:00
#PBS -q debug
#PBS -A GFS-DEV
#PBS -l select=4:ncpus=32
#PBS -V

set -x

# specify computation resource
export threads=1
export OMP_NUM_THREADS=$threads
export APRUN="mpiexec -l -n 128 -ppn 32 --cpu-bind core --depth 1"

echo "starting time"
date

############################################
# Loading module
############################################
module reset
module load intel/19.1.3.304
module load PrgEnv-intel/8.1.0
module load craype/2.7.8
module load cray-mpich/8.1.7
module load cray-pals/1.0.12
module load hdf5/1.10.6
module load netcdf/4.7.4
module load crtm/2.4.0.1
module load libjpeg/9c
module load prod_util/2.0.8
module load wgrib2/2.0.8
module list

msg="Starting 3drtma test"
postmsg "$logfile" "$msg"


export POSTGPEXEC=${svndir}/exec/upp.x

# specify your running and output directory
export startdate=2023040400
export fhr=000
export tmmark=tm00
export DATA=$rundir/3drtma_${startdate}

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
DateStr='${YY}-${MM}-${DD}_${HH}:00:00'
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
cp ${svndir}/parm/rrfs/postxconfig-NT-rrfs.txt postxconfig-NT.txt
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

${APRUN} ${POSTGPEXEC} < itag > wrfpost2.out

# operational 3drtma post processing generates 3 files
filelist="NATLEV00.tm00 \
          PRSLEV00.tm00"

for file in $filelist; do
export filein2=$file
ls -l ${filein2}
export err=$?

if [ $err = "0" ] ; then

 # use cmp to see if new pgb files are identical to the control one
 cmp ${filein2} $homedir/data_out/3drtma/${filein2}.${machine}

 # if not bit-identical, use cmp_grib2_grib2 to compare each grib record
 export err1=$?
 if [ $err1 -eq 0 ] ; then
  msg="3drtma test: your new post executable generates bit-identical ${filein2} as the develop branch"
  echo $msg
 else
  msg="3drtma test: your new post executable did not generate bit-identical ${filein2} as the develop branch"
  echo $msg
  $cmp_grib2_grib2 $homedir/data_out/3drtma/${filein2}.${machine} ${filein2} > ${filein2}.diff
 fi


else

 msg="3drtma test: post failed using your new post executable to generate ${filein2}"
 echo $msg 2>&1 | tee -a TEST_ERROR

fi
postmsg "$logfile" "$msg"
done

echo "PROGRAM IS COMPLETE!!!!!" 2>&1 | tee SUCCESS
msg="Ending 3drtma test"
postmsg "$logfile" "$msg"


