#!/bin/sh

#SBATCH -o out.post.fv3r_pe_test
#SBATCH -e out.post.fv3r_pe_test
#SBATCH -J fv3r_pe_test
#SBATCH -t 00:30:00
#SBATCH -N 5 --ntasks-per-node=12
#SBATCH -q batch
#SBATCH -A nems

set -x

# specify computation resource
export threads=1
export MP_LABELIO=yes
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

ulimit -s unlimited
export WGRIB2=wgrib2
export COMROOT=$rundir
#export CRTM_FIX=/apps/contrib/NCEPLIBS/orion/fix/crtm_v2.3.0

msg="Starting fv3r pe test"
postmsg "$logfile" "$msg"

export cmp_grib2_grib2=/home/wmeng/bin/cmp_grib2_grib2_new
# specify user's own post executable for testing
export POSTGPEXEC=${svndir}/exec/upp.x     


# specify forecast start time and hour for running your post job
export startdate=2023062800
export fhr=010

# specify your running and output directory
export DATA=$rundir/fv3r_${startdate}_pe_test
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
fileName='$homedir/data_in/fv3r/dynf${fhr}.nc'
IOFORM='netcdf'
grib='grib2'
DateStr='${YY}-${MM}-${DD}_${HH}:00:00'
MODELNAME='FV3R'
fileNameFlux='$homedir/data_in/fv3r/phyf${fhr}.nc'
/
&NAMPGB
KPO=47,PO=1000.,975.,950.,925.,900.,875.,850.,825.,800.,775.,750.,725.,700.,675.,650.,625.,600.,575.,550.,525.,500.,475.,450.,425.,400.,375.,350.,325.,300.,275.,250.,225.,200.,175.,150.,125.,100.,70.,50.,30.,20.,10.,7.,5.,3.,2.,1.,
/
EOF


rm -f fort.*

#cp /nwprod/nam.v3.1.16/fix/nam_micro_lookup.dat ./eta_micro_lookup.dat
cp ${svndir}/fix/nam_micro_lookup.dat ./eta_micro_lookup.dat

# copy flat files instead
cp ${svndir}/parm/postxconfig-NT-rrfs.txt ./postxconfig-NT.txt

cp ${svndir}/parm/params_grib2_tbl_new ./params_grib2_tbl_new

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

${APRUN} ${POSTGPEXEC} < itag > outpost_nems_${NEWDATE}

fhr=`expr $fhr + 0`
fhr2=`printf "%02d" $fhr`
#mv BGDAWP${fhr2}.tm00 BGDAWP${fhr2}.tm00.Grib2
#mv BGRD3D${fhr2}.tm00 BGRD3D${fhr2}.tm00.Grib2

filelist="PRSLEV${fhr2}.tm00 \
          NATLEV${fhr2}.tm00"

for file in $filelist; do
export filein2=$file
ls -l ${filein2}
export err=$?

if [ $err = "0" ] ; then

 # use cmp to see if new pgb files are identical to the control one
 cmp ${filein2} $homedir/data_out/fv3r/${filein2}.${machine}

 # if not bit-identical, use cmp_grib2_grib2 to compare each grib record
 export err1=$?
 if [ $err1 -eq 0 ] ; then
  msg="fv3r pe test: your new post executable generates bit-identical ${filein2} as the trunk"
  echo $msg
 else
  msg="fv3r pe test: your new post executable did not generate bit-identical ${filein2} as the trunk"
  echo $msg
  echo " start comparing each grib record and write the comparison result to *diff files"
  echo " check these *diff files to make sure your new post only change variables which you intend to change"
  $cmp_grib2_grib2 $homedir/data_out/fv3r/${filein2}.${machine} ${filein2} > ${filein2}.diff
 fi

else

 msg="fv3r pe test: post failed using your new post executable to generate ${filein2}"
 echo $msg 2>&1 | tee -a TEST_ERROR

fi
postmsg "$logfile" "$msg"
done

echo "PROGRAM IS COMPLETE!!!!!" 2>&1 | tee SUCCESS
msg="Ending fv3r test"
postmsg "$logfile" "$msg"
