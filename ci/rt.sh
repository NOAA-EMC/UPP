#!/bin/bash
######################################################################
# This script is desined for UPP regression tests run by UPP developer.
# Wen Meng, 12/2020, First version.
# Fernando Andrade-Maldonado 5/2023 rework for CLI Options
# Fernando Andrade-Maldonado / Wen Meng 9/2023 Add Hercules, fix typos, and refactor
#
######################################################################

git_branch="develop"
git_url="https://github.com/NOAA-EMC/UPP.git"
clone_on="no"

while getopts a:w:h:r:t:b:u:c opt; do
  case $opt in
    a) accnr=${OPTARG}
        ;;
    w) workdir=${OPTARG}
        ;;
    h) homedir=${OPTARG}
        ;;
    r) rundir=${OPTARG}
        ;;
    t) test_v=${OPTARG}
        ;;
    b) git_branch=${OPTARG}
        ;;
    u) git_url=${OPTARG}
        ;;
    c) clone_on="yes"
	;;
  esac
done

#UPP working copy
test_v=${test_v:-`pwd`/../}
if [[ $clone_on == "yes" ]]; then
  rm -rf $test_v
  mkdir -p $test_v
  git clone -b $git_branch $git_url $test_v
fi
export svndir=${test_v}

#Assume a nems account to run with
accnr=${accnr:-"nems"}

#Build UPP executable
build_exe=yes

#Choose run specific model
run_nmmb=yes
run_gfs=yes
run_gefs=yes
run_fv3r=yes
run_rap=yes
run_hrrr=yes
run_hafs=yes
run_rtma=yes

#find machine
mac=$(hostname | cut -c1-1)
mac2=$(hostname | cut -c1-2)
mac3=$(hostname | cut -c1-4)
if [ $mac2 = hf ]; then #for HERA
 export machine=HERA
 export homedir=${homedir:-"/scratch2/NAGAPE/epic/UPP/test_suite"}
 export rundir=${rundir:-"/scratch1/NCEPDEV/stmp2/${USER}"}
 module use /scratch1/NCEPDEV/nems/role.epic/hpc-stack/libs/intel-2022.1.2/modulefiles/stack
 module load hpc/1.2.0
 module load hpc-intel/2022.1.2
 module load hpc-impi/2022.1.2
 module load prod_util
elif [ $mac = O ] ; then
 export machine=ORION
 export homedir=${homedir:-"/work/noaa/epic/UPP"}
 export rundir=${rundir:-"/work2/noaa/stmp/$USER"}
 module use /work/noaa/epic-ps/role-epic-ps/hpc-stack/libs/intel-2022.1.2/modulefiles/stack
 module load hpc/1.2.0
 module load hpc-intel/2022.1.2
 module load hpc-impi/2022.1.2
 module load prod_util/1.2.2
elif [ $mac3 = herc ] ; then
 export machine=HERCULES
 export homedir=${homedir:-"/work/noaa/epic/UPP"}
 export rundir=${rundir:-"/work2/noaa/stmp/$USER"}
 module use /work/noaa/epic/role-epic/spack-stack/hercules/spack-stack-dev-20230717/envs/unified-env/install/modulefiles/Core
 module load stack-intel/2021.9.0
 module load stack-intel-oneapi-mpi/2021.9.0
 module load jasper/2.0.32
 module load prod-util/1.2.2
fi

#set working directory
export workdir=${workdir:-"`pwd`/work-upp-${machine}"}
rm -rf $workdir
mkdir -p $workdir

#differentiates for orion and hercules
export rundir="${rundir}/upp-${machine}"

#set log file
export logfile=`pwd`/rt.log.$machine
if [ -f $logfile ] ; then
 rm -r $logfile
fi
runtime_log=$homedir/scripts/runtime.log.$machine

#build executable
if [ "$build_exe" = "yes" ]; then
cd ${test_v}
mkdir -p ${test_v}/exec
cd ${test_v}/tests
./compile_upp.sh
status=$?
if [ $status -eq 0 ]; then
  msg="Building executable successfully"
else
  msg="Building executable with failure"
  exit
fi
postmsg "$logfile" "$msg"
fi

jobid_list=""

#execute nmmb grib2 test
if [ "$run_nmmb" = "yes" ]; then
cd $workdir
cp $homedir/jobs-dev/run_post_nmmb_Grib2_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr} run_post_nmmb_Grib2_${machine}.sh`
jobid_list=$jobid_list" "$job_id
cp $homedir/jobs-dev/run_post_nmmb_Grib2_pe_test_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr} run_post_nmmb_Grib2_pe_test_${machine}.sh`
jobid_list=$jobid_list" "${job_id}
fi

#execute fv3gefs test
if [ "$run_gefs" = "yes" ]; then
cd $workdir
cp $homedir/jobs-dev/run_post_fv3gefs_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr} run_post_fv3gefs_${machine}.sh`
jobid_list=$jobid_list" "${job_id}
cp $homedir/jobs-dev/run_post_fv3gefs_pe_test_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr} run_post_fv3gefs_pe_test_${machine}.sh`
jobid_list=$jobid_list" "${job_id}
fi

#execute rap test
if [ "$run_rap" = "yes" ]; then
cd $workdir
cp $homedir/jobs-dev/run_post_rap_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr} run_post_rap_${machine}.sh`
jobid_list=$jobid_list" "$job_id
cp $homedir/jobs-dev/run_post_rap_pe_test_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr} run_post_rap_pe_test_${machine}.sh`
jobid_list=$jobid_list" "${job_id}
fi

#execute hrrr test
if [ "$run_hrrr" = "yes" ]; then
cd $workdir
cp $homedir/jobs-dev/run_post_hrrr_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr} run_post_hrrr_${machine}.sh`
jobid_list=$jobid_list" "$job_id
cp $homedir/jobs-dev/run_post_hrrr_pe_test_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr} run_post_hrrr_pe_test_${machine}.sh`
jobid_list=$jobid_list" "${job_id}
fi

#execute fv3gfs test
if [ "$run_gfs" = "yes" ]; then
cd $workdir
cp $homedir/jobs-dev/run_post_fv3gfs_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr}  run_post_fv3gfs_${machine}.sh`
jobid_list=$jobid_list" "${job_id}
cp $homedir/jobs-dev/run_post_fv3gfs_pe_test_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr} run_post_fv3gfs_pe_test_${machine}.sh`
jobid_list=$jobid_list" "${job_id}
fi

#execute fv3r test
if [ "$run_fv3r" = "yes" ]; then
cd $workdir
cp $homedir/jobs-dev/run_post_fv3r_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr} run_post_fv3r_${machine}.sh`
jobid_list=$jobid_list" "${job_id}
cp $homedir/jobs-dev/run_post_fv3r_pe_test_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr} run_post_fv3r_pe_test_${machine}.sh`
jobid_list=$jobid_list" "${job_id}
fi

#execute fv3hafs test
if [ "$run_hafs" = "yes" ]; then
cd $workdir
cp $homedir/jobs-dev/run_post_fv3hafs_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr} run_post_fv3hafs_${machine}.sh`
jobid_list=$jobid_list" "${job_id}
cp $homedir/jobs-dev/run_post_fv3hafs_pe_test_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr} run_post_fv3hafs_pe_test_${machine}.sh`
jobid_list=$jobid_list" "${job_id}
fi

#execute rtma test
if [ "$run_rtma" = "yes" ]; then
cd $workdir
cp $homedir/jobs-dev/run_post_3drtma_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr} run_post_3drtma_${machine}.sh`
jobid_list=$jobid_list" "${job_id}
cp $homedir/jobs-dev/run_post_3drtma_pe_test_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr} run_post_3drtma_pe_test_${machine}.sh`
jobid_list=$jobid_list" "${job_id}
fi

echo "Job cards submitted for enabled tests, waiting on timestamps for finished jobs..."

#get run time for each test
sleep 30
for job_id in $jobid_list; do
  ic=1
  sleep_loop_max=300
  while [ $ic -le $sleep_loop_max ]; do
     job_id=`echo $job_id | cut -d"." -f1`
     status=`sacct --parsable -j $job_id --format=jobid,jobname,elapsed,state | cut -d"|" -f4|awk 'FNR == 2'`
     if [ "$status" = "COMPLETED" ]; then
       break
     else
      ic=`expr $ic + 1`
      sleep 15
     fi
  done
  if [ $ic -lt $sleep_loop_max ]; then
     runtime=`sacct --parsable -j $job_id --format=jobid,jobname,elapsed,state | cut -d"|" -f3|awk 'FNR == 2'`
     jobname=`sacct --parsable -j $job_id --format=jobid,jobname,elapsed,state | cut -d"|" -f2|awk 'FNR == 2'`
     runtime_b=`grep ${jobname} ${runtime_log} | awk '{print $2}' `
     echo "$runtime   $jobname ${runtime_b}"
     msg="Runtime: $jobname $runtime -- baseline ${runtime_b}"
     postmsg "$logfile" "$msg"
  fi
done

python ${test_v}/ci/rt-status.py

# Cleanup rt log
cd ${test_v}/ci
echo "rundir: ${rundir}" > rt.log.${machine}.temp
cat rt.log.${machine} | grep "test:" >> rt.log.${machine}.temp
cat rt.log.${machine} | grep "baseline" >> rt.log.${machine}.temp
cat rt.log.${machine}.temp > rt.log.${machine}
rm rt.log.${machine}.temp
mv rt.log.${machine} ${test_v}/tests/logs
