#!/bin/bash
######################################################################
# This script is desined for UPP regression tests run by UPP developer.
# Wen Meng, 12/2020, First version.
# Fernando Andrade-Maldonado 5/2023 rework for CLI Options
# Fernando Andrade-Maldonado / Wen Meng 9/2023 Add Hercules, fix typos, and refactor
# Fernando Andrade-Maldonado 4/2024 Additional Log info
######################################################################
set -xue
SECONDS=0

git_branch="develop"
git_url="https://github.com/NOAA-EMC/UPP.git"
clone_on="no"
disable_ifi="no" # don't use libIFI, even if it is present

while getopts a:w:h:r:t:b:u:cd opt; do
  case $opt in
    d) disable_ifi=yes
        ;;
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

if [[ -d $svndir/sorc/libIFI.fd/src/ ]] ; then
    have_ifi=yes
else
    have_ifi=no
fi

#Assume a nems account to run with
accnr=${accnr:-"rtrr"}

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

# Tests with IFI enabled only work if libIFI is present.
if [[ "$have_ifi" == yes && "$disable_ifi" == no ]] ; then
  run_hrrr_ifi=yes
  run_ifi_standalone_hrrr=yes
  run_fv3r_ifi=yes
  run_ifi_standalone_fv3r=yes
else
  # Cannot run these without ifi
  run_hrrr_ifi=no
  run_ifi_standalone_hrrr=no
  run_fv3r_ifi=no
  run_ifi_standalone_fv3r=no
fi

#find machine
mac=$(hostname | cut -c1-1)
mac2=$(hostname | cut -c1-2)
mac3=$(hostname | cut -c1-4)
if [ $mac2 = hf ]; then # for HERA
 export machine=HERA
 export homedir=${homedir:-"/scratch2/NAGAPE/epic/UPP/test_suite"}
 export rundir=${rundir:-"/scratch1/NCEPDEV/stmp2/${USER}"}
 module use /contrib/spack-stack/spack-stack-1.8.0/envs/ue-intel-2021.5.0/install/modulefiles/Core
 module load stack-intel/2021.5.0
 module load stack-intel-oneapi-mpi/2021.5.1
 module load prod_util/2.1.1
elif [ $mac3 = orio ] ; then
 export machine=ORION
 export homedir=${homedir:-"/work/noaa/epic/UPP"}
 export rundir=${rundir:-"/work2/noaa/stmp/$USER"}
 module use /apps/contrib/spack-stack/spack-stack-1.8.0/envs/ue-intel-2021.9.0/install/modulefiles/Core
 module load stack-intel/2021.9.0
 module load stack-intel-oneapi-mpi/2021.9.0
 module load prod_util/2.1.1
 module load python/3.10.8
elif [ $mac3 = herc ] ; then
 export machine=HERCULES
 export homedir=${homedir:-"/work/noaa/epic/UPP"}
 export rundir=${rundir:-"/work2/noaa/stmp/$USER"}
 module use /apps/contrib/spack-stack/spack-stack-1.8.0/envs/ue-intel-2021.9.0/install/modulefiles/Core
 module load stack-intel/2021.9.0
 module load stack-intel-oneapi-mpi/2021.9.0
 module load prod_util/2.1.1
 module load python/3.10.8
fi

#set working directory
export workdir=${workdir:-"`pwd`/work-upp-${machine}"}
rm -rf $workdir
mkdir -p $workdir

#differentiates for orion and hercules
export rundir="${rundir}/upp-${machine}"
test -d "${rundir}" || mkdir -p "${rundir}"

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
  ./compile_upp.sh -o upp_no_ifi.x
  status=$?
  if [ $status -eq 0 ]; then
    msg="Building executable successfully"
  else
    msg="Building executable with failure"
    postmsg "$logfile" "$msg"
    exit 2
  fi

  if [[ "$have_ifi" == yes && "$disable_ifi" == no ]] ; then
    ./compile_upp.sh -a -o upp_with_ifi.x -I -B
    status=$?
    if [ $status -eq 0 ]; then
      msg="Building UPP+IFI executables successfully"
    else
      msg="Building UPP+IFI executables with failure"
      postmsg "$logfile" "$msg"
      exit 2
    fi
    ln -s upp_with_ifi.x $svndir/exec/upp.x
  else
    ln -s upp_no_ifi.x $svndir/exec/upp.x
  fi

  postmsg "$logfile" "$msg"
fi

jobid_list=""
set -xe
#execute ifi tests           
if [ "${run_hrrr_ifi:-no}" = "yes" ]; then
cd $workdir
cp $svndir/ci/jobs-dev/run_post_hrrr_ifi_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr} run_post_hrrr_ifi_${machine}.sh`
jobid_list=$jobid_list" "${job_id}
dep_job_id=$job_id
  if [ "$run_ifi_standalone_hrrr" = "yes" ]; then
    cp $svndir/ci/jobs-dev/run_ifi_standalone_hrrr_${machine}.sh .
    job_id=`sbatch --parsable -A ${accnr} --dependency=afterany:$dep_job_id run_ifi_standalone_hrrr_${machine}.sh`
    jobid_list=$jobid_list" "${job_id}
  fi
fi

if [ "$run_fv3r_ifi" = "yes" ]; then
cd $workdir
cp $svndir/ci/jobs-dev/run_post_fv3r_ifi_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr} run_post_fv3r_ifi_${machine}.sh`
jobid_list=$jobid_list" "${job_id}
dep_job_id=$job_id
  if [ "$run_ifi_standalone_fv3r" = "yes" ]; then
    cp $svndir/ci/jobs-dev/run_ifi_standalone_fv3r_${machine}.sh .
    job_id=`sbatch --parsable -A ${accnr} --dependency=afterany:$dep_job_id run_ifi_standalone_fv3r_${machine}.sh`
    jobid_list=$jobid_list" "${job_id}
  fi
fi

#execute nmmb grib2 test
if [ "$run_nmmb" = "yes" ]; then
cd $workdir
cp $svndir/ci/jobs-dev/run_post_nmmb_Grib2_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr} run_post_nmmb_Grib2_${machine}.sh`
jobid_list=$jobid_list" "$job_id
cp $svndir/ci/jobs-dev/run_post_nmmb_Grib2_pe_test_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr} run_post_nmmb_Grib2_pe_test_${machine}.sh`
jobid_list=$jobid_list" "${job_id}
fi

#execute fv3gefs test
if [ "$run_gefs" = "yes" ]; then
cd $workdir
cp $svndir/ci/jobs-dev/run_post_fv3gefs_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr} run_post_fv3gefs_${machine}.sh`
jobid_list=$jobid_list" "${job_id}
cp $svndir/ci/jobs-dev/run_post_fv3gefs_pe_test_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr} run_post_fv3gefs_pe_test_${machine}.sh`
jobid_list=$jobid_list" "${job_id}
fi

#execute rap test
if [ "$run_rap" = "yes" ]; then
cd $workdir
cp $svndir/ci/jobs-dev/run_post_rap_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr} run_post_rap_${machine}.sh`
jobid_list=$jobid_list" "$job_id
cp $svndir/ci/jobs-dev/run_post_rap_pe_test_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr} run_post_rap_pe_test_${machine}.sh`
jobid_list=$jobid_list" "${job_id}
fi

#execute hrrr test
if [ "$run_hrrr" = "yes" ]; then
cd $workdir
cp $svndir/ci/jobs-dev/run_post_hrrr_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr} run_post_hrrr_${machine}.sh`
jobid_list=$jobid_list" "$job_id
cp $svndir/ci/jobs-dev/run_post_hrrr_pe_test_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr} run_post_hrrr_pe_test_${machine}.sh`
jobid_list=$jobid_list" "${job_id}
fi

#execute fv3gfs test
if [ "$run_gfs" = "yes" ]; then
cd $workdir
cp $svndir/ci/jobs-dev/run_post_fv3gfs_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr}  run_post_fv3gfs_${machine}.sh`
jobid_list=$jobid_list" "${job_id}
cp $svndir/ci/jobs-dev/run_post_fv3gfs_pe_test_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr} run_post_fv3gfs_pe_test_${machine}.sh`
jobid_list=$jobid_list" "${job_id}
fi

#execute fv3r test
if [ "$run_fv3r" = "yes" ]; then
cd $workdir
cp $svndir/ci/jobs-dev/run_post_fv3r_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr} run_post_fv3r_${machine}.sh`
jobid_list=$jobid_list" "${job_id}
cp $svndir/ci/jobs-dev/run_post_fv3r_pe_test_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr} run_post_fv3r_pe_test_${machine}.sh`
jobid_list=$jobid_list" "${job_id}
cp $svndir/ci/jobs-dev/run_post_fv3r_ifi_missing_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr} run_post_fv3r_ifi_missing_${machine}.sh`
jobid_list=$jobid_list" "${job_id}
fi

#execute fv3hafs test
if [ "$run_hafs" = "yes" ]; then
cd $workdir
cp $svndir/ci/jobs-dev/run_post_fv3hafs_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr} run_post_fv3hafs_${machine}.sh`
jobid_list=$jobid_list" "${job_id}
cp $svndir/ci/jobs-dev/run_post_fv3hafs_pe_test_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr} run_post_fv3hafs_pe_test_${machine}.sh`
jobid_list=$jobid_list" "${job_id}
fi

#execute rtma test
if [ "$run_rtma" = "yes" ]; then
cd $workdir
cp $svndir/ci/jobs-dev/run_post_3drtma_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr} run_post_3drtma_${machine}.sh`
jobid_list=$jobid_list" "${job_id}
cp $svndir/ci/jobs-dev/run_post_3drtma_pe_test_${machine}.sh .
job_id=`sbatch --parsable -A ${accnr} run_post_3drtma_pe_test_${machine}.sh`
jobid_list=$jobid_list" "${job_id}
fi
set +xe
echo "Job cards submitted for enabled tests, waiting on timestamps for finished jobs..."

#get run time for each test
some_failed=NO
sleep 30
for job_id in $jobid_list; do
  ic=1
  sleep_loop_max=300
  while [ $ic -le $sleep_loop_max ]; do
     job_id=`echo $job_id | cut -d"." -f1`
     status=`sacct --parsable -j $job_id --format=jobid,jobname,elapsed,state | cut -d"|" -f4|awk 'FNR == 2'`
     if [ "$status" = "COMPLETED" ]; then
       break
     elif ( echo "$status" | grep -E 'FAIL|TIMEOUT|CANCEL|DEAD|SIGNAL|SPECIAL' > /dev/null ) ; then
       some_failed=YES
       break
     else
      ic=`expr $ic + 1`
      sleep 15
     fi
  done
  if [ $ic -lt $sleep_loop_max ]; then
     runtime=`sacct --parsable -j $job_id --format=jobid,jobname,elapsed,state | cut -d"|" -f3|awk 'FNR == 2'`
     jobname=`sacct --parsable -j $job_id --format=jobid,jobname,elapsed,state | cut -d"|" -f2|awk 'FNR == 2'`
     runtime_b=`grep "^${jobname}" ${runtime_log} | awk '{print $2}'`
     echo "$runtime   $jobname ${runtime_b}"
     msg="Runtime: $jobname $runtime -- baseline ${runtime_b}"
     postmsg "$logfile" "$msg"
  fi
done

elapsed_time=$( printf '%02dh:%02dm:%02ds\n' $((SECONDS%86400/3600)) $((SECONDS%3600/60)) $((SECONDS%60)) )

python ${test_v}/ci/rt-status.py
test_results=$?

if [ $some_failed = YES ] ; then
  test_results=99
  echo WARNING: some tests exited with non-zero status.
fi

# Cleanup rt log
cd ${test_v}

UPP_HASH=$(git rev-parse HEAD)
SUBMODULE_HASHES=$(git submodule status --recursive)
DATE="$(date '+%Y%m%d %T')"

cd ${test_v}/ci

cat << EOF > rt.log.${machine}.temp
===== Start of UPP Regression Testing Log =====
UPP Hash Tested:
${UPP_HASH}

Submodule hashes:
${SUBMODULE_HASHES}

Run directory: ${rundir}
Baseline directory: ${homedir}

Total runtime: ${elapsed_time}
Test Date: ${DATE}
Summary Results:

EOF


if [ $some_failed = YES ] ; then
    echo "Warning: some tests exited with non-zero. status" >> rt.log.${machine}.temp
    echo >> rt.log.${machine}.temp
fi

cat rt.log.${machine} | grep "test:" >> rt.log.${machine}.temp
cat rt.log.${machine} | grep "baseline" >> rt.log.${machine}.temp
python ${test_v}/ci/rt-status.py >> rt.log.${machine}.temp
echo "===== End of UPP Regression Testing Log =====" >> rt.log.${machine}.temp
mv rt.log.${machine}.temp rt.log.${machine}
mv rt.log.${machine} ${test_v}/tests/logs

# should indicate failure to Jenkins
if [ $test_results -ne 0 ]; then
   python ${test_v}/ci/rt-status.py > changed_results.txt
   if [ $some_failed = YES ]; then
     echo "Warning: some tests exited with non-zero status." >> changed_results.txt
   fi
   exit 1
fi
