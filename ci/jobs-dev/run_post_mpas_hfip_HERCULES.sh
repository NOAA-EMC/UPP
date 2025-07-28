#!/bin/bash

#SBATCH -o out.post.mpas_hfip
#SBATCH -e out.post.mpas_hfip
#SBATCH -J mpas_hfip_test 
#SBATCH -t 00:30:00
#SBATCH --ntasks=256
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH -q batch
#SBATCH -A nems
#SBATCH --exclusive

set -x

# specify computation resources
export threads=4
export MP_LABELIO=yes
export OMP_NUM_THREADS=$threads
export APRUN="srun"

echo "starting time"
date

############################################
# Loading modules
############################################
module use ${svndir}/modulefiles
module load hercules_$compiler
module load prod_util/2.1.1
module load wgrib2/3.6.0
module list

ulimit -s unlimited

msg="Starting mpas_hfip test"
postmsg "$logfile" "$msg"

export POSTGPEXEC=${svndir}/exec/upp.x

# specify forecast start time and hour
export startdate=2024100700
export fhr=048

# specify your running and output directory
export DATA=$rundir/mpas_hfip_${startdate}
rm -rf $DATA; mkdir -p $DATA
cd $DATA

export NEWDATE=`${NDATE} +${fhr} $startdate`
export YY=`echo $NEWDATE | cut -c1-4`
export MM=`echo $NEWDATE | cut -c5-6`
export DD=`echo $NEWDATE | cut -c7-8`
export HH=`echo $NEWDATE | cut -c9-10`

cat > itag <<EOF
&model_inputs
    fileName='$homedir/data_in/mpas_hfip/MPAS-A_out.${YY}-${MM}-${DD}_${HH}.00.00.nc'
    ioform = 'netcdfpara'
    grib = 'grib2'
    datestr = '${YY}-${MM}-${DD}_${HH}:00:00'
    modelname = 'RAPR'
    submodelname = 'MPAS'
/
EOF

# copy fix data
cp ${svndir}/fix/rap_micro_lookup.dat .
cp ${svndir}/fix/nam_micro_lookup.dat .
cp ${svndir}/parm/mpas/postxconfig-NT-hfip_mpas.txt ./postxconfig-NT.txt
cp ${svndir}/parm/params_grib2_tbl_new ./params_grib2_tbl_new

# get crtm fix files
for what in \
    "FASTEM4.MWwater" "FASTEM5.MWwater" "FASTEM6.MWwater" "NPOESS.IRice" "NPOESS.IRland" \
    "NPOESS.IRsnow" "Nalli.IRwater" "abi_gr" "ahi_himawari8" "amsre_aqua" \
    "imgr_g11" "imgr_g12" "imgr_g13" "imgr_g15" "imgr_insat3d" "imgr_mt1r" "imgr_mt2" \
    "ssmi_f13" "ssmi_f14" "ssmi_f15" "ssmis_f16" "ssmis_f17" "ssmis_f18" "ssmis_f19" "ssmis_f20" \
    "seviri_m10" "tmi_trmm" "v.seviri_m10" ; do
  for coef in Spc Tau ; do
    file="${CRTM_FIX}/${what}.${coef}Coeff.bin"
    if [[ -s "$file" ]] ; then
      ln -s "$file" .
    fi
  done
done

for what in 'Aerosol' 'Cloud' ; do
    ln -s "${CRTM_FIX}/${what}Coeff.bin" .
done

for what in  ${CRTM_FIX}/*Emis* ; do
   ln -s $what .
done

# Run the UPP
export PGBOUT=pgbfile
${APRUN} ${POSTGPEXEC} < itag > outpost_mpas_hfip_${startdate}

################################################
# Compare with baseline data
################################################
fhr=$((10#$fhr))
fhr2=$(printf "%02d" "$fhr")

# MPAS_HFIP post processing generates 3 files
filelist="NATLEV.GrbF$fhr2 \
          PRSLEV.GrbF$fhr2 \
          2DFLD.GrbF$fhr2"

for file in $filelist; do
export filein2=$file
ls -l ${filein2}
export err=$?

if [ $err = "0" ] ; then
 # use cmp to see if new pgb files are identical to the control one
 cmp ${filein2} $homedir/data_out_$compiler/mpas_hfip/${filein2}.${machine}

 # if not bit-identical, use cmp_grib2_grib2 to compare each grib record
 export err1=$?
 if [ $err1 -eq 0 ] ; then
  msg="mpas_hfip test: your new post executable generates bit-identical ${filein2} as the develop branch"
  echo $msg
 else
  msg="mpas_hfip test: your new post executable did not generate bit-identical ${filein2} as the develop branch"
  echo $msg
  echo " start comparing each grib record and write the comparison result to *diff files"
  echo " check these *diff files to make sure your new post only change variables which you intend to change"
  $cmp_grib2_grib2 $homedir/data_out_$compiler/mpas_hfip/${filein2}.${machine} ${filein2} > ${filein2}.diff
 fi
else
 msg="mpas_hfip test: post failed using your new post executable to generate ${filein2}"
 echo $msg 2>&1 | tee -a TEST_ERROR
fi

postmsg "$logfile" "$msg"
done

echo "PROGRAM IS COMPLETE!!!!!" 2>&1 | tee SUCCESS
msg="Ending mpas_hfip test"
postmsg "$logfile" "$msg"
