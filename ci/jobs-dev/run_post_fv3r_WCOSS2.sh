#!/bin/bash

#PBS -o out.post.fv3r
#PBS -e out.post.fv3r
#PBS -N fv3r.test
#PBS -l walltime=00:30:00
#PBS -q debug
#PBS -A GFS-DEV
#PBS -l place=vscatter,select=2:ncpus=48
#PBS -V

set -x

# specify computation resource
export threads=1
export OMP_NUM_THREADS=$threads
export APRUN="mpiexec -l -n 96 -ppn 48"

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
module load libjpeg/9c
module load prod_util/2.0.8
module load crtm/2.4.0.1
module list

msg="Starting fv3r test"
postmsg "$logfile" "$msg"

export cmp_grib2_grib2=/u/wen.meng/bin/cmp_grib2_grib2_new
export POSTGPEXEC=${svndir}/exec/upp.x

# specify forecast start time and hour for running your post job
export startdate=2023062800
export fhr=010

# specify your running and output directory
export DATA=$rundir/fv3r_${startdate}
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

#copy fix data
cp ${svndir}/fix/nam_micro_lookup.dat ./eta_micro_lookup.dat
cp ${svndir}/parm/rrfs/postxconfig-NT-rrfs.txt ./postxconfig-NT.txt
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
  msg="fv3r test: your new post executable generates bit-identical ${filein2} as the develop branch"
  echo $msg
 else
  msg="fv3r test: your new post executable did not generate bit-identical ${filein2} as the develop branch"
  echo $msg
  echo " start comparing each grib record and write the comparison result to *diff files"
  echo " check these *diff files to make sure your new post only change variables which you intend to change"
  $cmp_grib2_grib2 $homedir/data_out/fv3r/${filein2}.${machine} ${filein2} > ${filein2}.diff
 fi

else

 msg="fv3r test: post failed using your new post executable to generate ${filein2}"
 echo $msg 2>&1 | tee -a TEST_ERROR

fi
postmsg "$logfile" "$msg"
done

echo "PROGRAM IS COMPLETE!!!!!" 2>&1 | tee SUCCESS
msg="Ending fv3r test"
postmsg "$logfile" "$msg"
