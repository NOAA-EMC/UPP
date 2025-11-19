#!/bin/bash

#PBS -o out.post.gfs
#PBS -e out.post.gfs
#PBS -N gfs_test
#PBS -l walltime=00:40:00
#PBS -q dev
#PBS -A GFS-DEV
#PBS -l place=vscatter,select=4:ncpus=48:mem=300GB
#PBS -V

set -x

# specify computation resources
export threads=1
export OMP_NUM_THREADS=$threads
export APRUN="mpiexec -l -n 192 -ppn 48"

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

msg="Starting gfs test"
postmsg "$logfile" "$msg"

export POSTGPEXEC=$svndir/exec/upp.x

# specify forecast start time and hour
export startdate=2025012600
export fhr=006
export cyc=`echo $startdate |cut -c9-10`

# specify your running and output directory
export DATA=$rundir/gfs_${startdate}
rm -rf $DATA; mkdir -p $DATA
cd $DATA

export NEWDATE=`${NDATE} +${fhr} $startdate`
export YY=`echo $NEWDATE | cut -c1-4`
export MM=`echo $NEWDATE | cut -c5-6`
export DD=`echo $NEWDATE | cut -c7-8`
export HH=`echo $NEWDATE | cut -c9-10`

cat > itag <<EOF
&model_inputs
fileName='$homedir/data_in/gfs/gfs.t${cyc}z.atmf${fhr}.nc'
IOFORM='netcdf'
grib='grib2'
DateStr='${YY}-${MM}-${DD}_${HH}:00:00'
MODELNAME='GFS'
fileNameFlux='$homedir/data_in/gfs/gfs.t${cyc}z.sfcf${fhr}.nc'
/
&NAMPGB
KPO=57,PO=1000.,975.,950.,925.,900.,875.,850.,825.,800.,775.,750.,725.,700.,675.,650.,625.,600.,575.,550.,525.,500.,475.,450.,425.,400.,375.,350.,325.,300.,275.,250.,225.,200.,175.,150.,125.,100.,70.,50.,40.,30.,20.,15.,10.,7.,5.,3.,2.,1.,0.7,0.4,0.2,0.1,0.07,0.04,0.02,0.01,rdaod=.true.,
/
EOF

# copy fix data
cp ${svndir}/fix/nam_micro_lookup.dat ./eta_micro_lookup.dat
cp ${svndir}/parm/params_grib2_tbl_new ./params_grib2_tbl_new

# get crtm fix files
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

# Generate master and flux files
cp ${svndir}/parm/gfs/postxconfig-NT-gfs-two.txt ./postxconfig-NT.txt
${APRUN} ${POSTGPEXEC} < itag > outpost_gfs_master_${NEWDATE}

# Generate goes file
cp ${svndir}/parm/gfs/postxconfig-NT-gfs-goes.txt ./postxconfig-NT.txt
${APRUN} ${POSTGPEXEC} < itag > outpost_gfs_goes_${NEWDATE}

################################################
# Compare with baseline data
################################################
fhr=$((10#$fhr))
FH3=$(printf "%03d" "$fhr")
FH2=$(printf "%02d" "$fhr")
mv GFSPRS.GrbF${FH2} gfs.t${cyc}z.master.grb2f${FH3}
mv GFSFLX.GrbF${FH2} gfs.t${cyc}z.sfluxgrbf${FH3}.grib2
mv GFSGOES.GrbF${FH2} gfs.t${cyc}z.special.grb2f${FH3}

# GFS post processing generates 3 files
filelist="gfs.t${cyc}z.master.grb2f${FH3} \
          gfs.t${cyc}z.sfluxgrbf${FH3}.grib2 \
          gfs.t${cyc}z.special.grb2f${FH3} "

for file in $filelist; do

export filein2=$file
ls -l ${filein2}
export err=$?

if [ $err = "0" ] ; then
 # use cmp to see if new pgb files are identical to the control one
 cmp ${filein2} $homedir/data_out/gfs/${filein2}.${machine}

 # if not bit-identical, use cmp_grib2_grib2 to compare each grib record

 export err1=$?
 if [ $err1 -eq 0 ] ; then
  msg="gfs test: your new post executable generates bit-identical ${filein2} as the develop branch"
  echo $msg
 else
  msg="gfs test: your new post executable did not generate bit-identical ${filein2} as the develop branch"
  echo $msg
  echo " start comparing each grib record and write the comparison result to *diff files"
  echo " check these *diff files to make sure your new post only change variables which you intend to change"
  $cmp_grib2_grib2 $homedir/data_out/gfs/${filein2}.${machine} ${filein2} > ${filein2}.diff
 fi
else
 msg="gfs test: post failed using your new post executable to generate ${filein2}"
 echo $msg 2>&1 | tee -a TEST_ERROR
fi

postmsg "$logfile" "$msg"
done

echo $?
echo "PROGRAM IS COMPLETE!!!!!!" 2>&1 | tee SUCCESS
msg="Ending gfs test"
postmsg "$logfile" "$msg"
