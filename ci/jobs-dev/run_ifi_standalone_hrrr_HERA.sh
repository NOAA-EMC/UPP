#!/bin/sh

#SBATCH -o out.post.ifi_standalone_hrrr
#SBATCH -e out.post.ifi_standalone_hrrr
#SBATCH -J ifi_standalone_hrrr_test
#SBATCH -t 00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --exclusive
#SBATCH --partition bigmem
#SBATCH -q batch
#SBATCH -A ovp

# specify computation resource
export threads=40
export MP_LABELIO=yes
export OMP_NUM_THREADS=$threads
export OMP_STACKSIZE=512M
export APRUN="srun"

############################################
# Loading module
############################################

# EXPORT list here

module purge
module use /contrib/spack-stack/spack-stack-1.8.0/envs/ue-intel-2021.5.0/install/modulefiles/Core
module load stack-intel/2021.5.0
module load stack-intel-oneapi-mpi/2021.5.1
module load libpng/1.6.37
module load jasper/2.0.32
module load prod_util/2.1.1
module load crtm/2.4.0.1
module load nccmp
module load netcdf-cxx4/4.3.1
module list

ulimit -s unlimited
ulimit

msg="Starting ifi_standalone_hrrr test"
postmsg "$logfile" "$msg"

export cmp_grib2_grib2=/home/Wen.Meng/bin/cmp_grib2_grib2_new
FIPEXEC=${svndir}/exec/fip2-lookalike.x

# use the UPP run directory so we get the input files in the expected format
export startdate=2020060118
export DATA=$rundir/hrrr_ifi_${startdate}
cd $DATA

upp_output=cat_vars_0.nc
ifi_standalone_output=icing-category-output.nc
diff_file=cat_vars_0.nc.diff

$APRUN --cpus-per-task=$OMP_NUM_THREADS --nodes=1 --ntasks=1 --exclusive \
     "$FIPEXEC" -u hybr_vars_0.nc hybr_vars_0.nc .

nccmp -n 20 -dfc1 -v ICE_PROB,ICE_SEV_CAT,SLD,WMO_ICE_SEV_CAT "$upp_output" "$ifi_standalone_output" 2>&1 | tee "$diff_file"
export err1=$?

if [ -s "$ifi_standalone_output" ] ; then
 if [ $err1 -eq 0 ] && ! [ -s "$diff_file" ] ; then
   msg="ifi standalone_fv3r test: Passed. libIFI standalone program and IFI in UPP produce identical results"
   rm -f "$diff_file"
   echo $msg
 else
   msg="ifi standalone_fv3r test: Failed. Differences detected between libIFI standalone and UPP IFI output."
   echo $msg
 fi
else
  msg="ifi standalone_fv3r test: Failed. ifi standalone failed using your new executable to generate $ifi_standalone_output"
  echo $msg 2>&1 | tee -a TEST_ERROR
fi
postmsg "$logfile" "$msg"

echo "PROGRAM IS COMPLETE!!!!!" 2>&1 | tee SUCCESS
msg="Ending ifi_standalone_hrrr test"
postmsg "$logfile" "$msg"
