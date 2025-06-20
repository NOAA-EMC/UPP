#!/bin/bash

#PBS -o out.post.fv3gfs
#PBS -e out.post.fv3gfs
#PBS -N fv3gfs.test
#PBS -l walltime=00:30:00
#PBS -q debug
#PBS -A GFS-DEV
#PBS -l place=vscatter,select=2:ncpus=48:mem=300GB
#PBS -V

set -x

# specify computation resource
export threads=1
export OMP_NUM_THREADS=$threads
export APRUN="mpiexec -l -n 96 -ppn 48"
export APRUN_DWN="mpiexec -l -np 24 --cpu-bind verbose,core cfp"

############################################
# Loading module
############################################
module reset
module load intel/19.1.3.304
module load PrgEnv-intel/8.1.0
module load craype/2.7.8
module load cray-mpich/8.1.7
module load cray-pals/1.0.12
module load cfp/2.0.4
module load hdf5/1.10.6
module load netcdf/4.7.4
module load crtm/2.4.0.1
module load prod_util/2.0.10
module load prod_envir/2.0.8
module load libjpeg/9c
module load grib_util/1.2.2
module load wgrib2/2.0.8
module list

#export WGRIB2=wgrib2

msg="Starting fv3gfs test"
postmsg "$logfile" "$msg"


export POSTGPEXEC=$svndir/exec/upp.x

# specify forecast start time and hour
export startdate=2019083000
export fhr=006
export cyc=`echo $startdate |cut -c9-10`

# specify your running and output directory
export DATA=$rundir/fv3gfs_${startdate}
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

#copy fix data
cp ${svndir}/fix/nam_micro_lookup.dat ./eta_micro_lookup.dat
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

#Generate master and flux files
cp ${svndir}/parm/gfs/postxconfig-NT-gfs-two.txt ./postxconfig-NT.txt
${APRUN} ${POSTGPEXEC} < itag > outpost_master_${NEWDATE}

#Generate goes file
cp ${svndir}/parm/gfs/postxconfig-NT-gfs-goes.txt ./postxconfig-NT.txt
${APRUN} ${POSTGPEXEC} < itag > outpost_goes_${NEWDATE}

#############################################################

################################################
# Compare with baseline data
################################################
FH3=$(printf %03i $fhr)
FH2=$(printf %02i $fhr)
mv GFSPRS.GrbF${FH2} gfs.t${cyc}z.master.grb2f${FH3}
mv GFSFLX.GrbF${FH2} gfs.t${cyc}z.sfluxgrbf${FH3}.grib2
mv GFSGOES.GrbF${FH2} gfs.t${cyc}z.special.grb2f${FH3}

# compare master, flux, simulated goes and wafs files
filelist="gfs.t${cyc}z.master.grb2f${FH3} \
          gfs.t${cyc}z.sfluxgrbf${FH3}.grib2 \
          gfs.t${cyc}z.special.grb2f${FH3} "
#          gfs.t${cyc}z.wafs.grb2f${post_times} "

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
  msg="fv3gfs test: your new post executable generates bit-identical ${filein2} as the develop branch"
  echo $msg
 else
  msg="fv3gfs test: your new post executable did not generate bit-identical ${filein2} as the develop branch"
  echo $msg
  echo " start comparing each grib record and write the comparison result to *diff files"
  echo " check these *diff files to make sure your new post only change variables which you intend to change"

  # compare grib message via cmp_grib2_grib2
  $cmp_grib2_grib2 $homedir/data_out/gfs/${filein2}.${machine} ${filein2} > ${filein2}.diff
 fi

else

 msg="fv3gfs test: post failed using your new post executable to generate ${filein2}"
 echo $msg 2>&1 | tee -a TEST_ERROR

fi
postmsg "$logfile" "$msg"
done

echo $?
echo "PROGRAM IS COMPLETE!!!!!!" 2>&1 | tee SUCCESS
msg="Ending fv3gfs test"
postmsg "$logfile" "$msg"
