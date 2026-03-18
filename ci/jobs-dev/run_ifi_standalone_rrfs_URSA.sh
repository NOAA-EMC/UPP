#!/bin/bash

#SBATCH -o out.post.ifi_standalone_rrfs
#SBATCH -e out.post.ifi_standalone_rrfs
#SBATCH -J ifi_standalone_rrfs_test
#SBATCH -t 00:30:00
#SBATCH --ntasks 240
#SBATCH --tasks-per-node 48
#SBATCH -q batch
#SBATCH -A rtrr
#SBATCH --exclusive

set -x

# specify computation resources
export threads=40
export MP_LABELIO=yes
export OMP_NUM_THREADS=$threads
export APRUN="srun"

echo "starting time"
date

############################################
# Loading modules
############################################
module purge
module use $svndir/modulefiles
module load ursa_$compiler
module load ursa_${compiler}_ifi_test_prereqs
module load wgrib2/3.6.0
module load prod_util/2.1.1
module load nccmp/1.9.1.0
module list

ulimit -s unlimited
ulimit

msg="Starting ifi_standalone_rrfs test"
postmsg "$logfile" "$msg"


FIPEXEC=${svndir}/exec/fip_runner

# use the UPP run directory so we get the input files in the expected format
export startdate=2025040112
export DATA=$rundir/rrfs_ifi_${startdate}
cd $DATA

upp_output=cat_vars_0.nc
ifi_standalone_output=20250401/fip_icing_category.20250401_g_120000_f_00064800.nc
diff_file=cat_vars_0.nc.diff

$APRUN --cpus-per-task=$OMP_NUM_THREADS --nodes=1 --ntasks=1 --exclusive \
     "$FIPEXEC" -u hybr_vars_0.nc hybr_vars_0.nc .

nccmp -dfc1 -v ICE_PROB,ICE_SEV_CAT,SLD,WMO_ICE_SEV_CAT "$upp_output" "$ifi_standalone_output" 2>&1 | tee "$diff_file"
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
msg="Ending ifi_standalone_rrfs test"
postmsg "$logfile" "$msg"
