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
print_full_help="no"

usage() {
  set +xue

  if [[ "$#" -gt 0 ]] ; then
    echo
    echo "------------------------------------------------------------------------"
    echo "$@"
    echo "------------------------------------------------------------------------"
  fi

  cat<<EOF

Synopsis: rt.sh -a account -r /path/to/scrub/space [-options] [compiler]
Executes UPP regression tests. Includes IFI tests if ../sorc/libIFI.fd exists.

Results are here:
  ../tests/logs/MACHINE_compiler.log = report of regression tests for each machine and compiler.
  changed_results.txt = A list of tests whose results have changed.

Always set these:
  -a account = accounting code for job submission. Default account is often overused. Always set this!
  -r rundir = path to a scrub space. Default area is often over quota. Always set this!

General options:
  -d = disable ifi tests even if ifi is available
  -h homedir = path to the regression test data
  -w workdir = directory with per-job batch and log files.
  -H = print full help message including special-use option flags.
EOF

  if [[ "$print_full_help" == YES ]] ; then
cat<<EOF

Special run mode: run rt.sh outside the repository. Automatically clones the repository.
Syntax: rt.sh -a account -r /path/to/scrub/space -c -u url -b branch [options] [compiler]

Additional options:
  -c = Tells rt.sh it is running outside a repository.
  -t test_v = Location to clone the repository. Default: Overwrite .. with the clone.
  -u url = Mandatory: URL of a repository to clone. Not for general use.
  -b branch = Mandatory: branch in the repository to clone
EOF
  fi
}

set +x
export OPTERR=1
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
    H) print_full_help=YES ; usage ; exit 1
        ;;
    *)
       usage FATAL ERROR: Invalid -option. See error message above. 1>&2
       exit 2
        ;;
  esac
done
set -x

if [[ "$OPTIND" > "$#" ]] ; then
  compiler=MISSING
else
  shift $(( OPTIND - 1 ))
  if [[ "$#" -gt 1 ]] ; then
    echo "ERROR: Expected at most 1 positional argument but found:" "$@" 1>&2
    usage FATAL ERROR: too many arguments. See error message above. 1>&2
    exit 2
  fi
  compiler="$1"
fi

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
elif [ $mac2 = uf ]; then # for Ursa
 export machine=URSA
 export homedir=${homedir:-"/scratch3/BMC/wrfruc/Samuel.Trahan/upp-ursa/test_suite"}
 export rundir=${rundir:-"/scratch3/BMC/wrfruc/Samuel.Trahan/scrub"}
 module use /contrib/spack-stack/spack-stack-1.9.1/envs/ue-oneapi-2024.2.1/install/modulefiles/Core
 module load stack-oneapi/2024.2.1
 module load stack-intel-oneapi-mpi/2021.13
 module load prod_util/2.1.1
 module load python/3.11.7
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

if [[ "$compiler" == MISSING ]] ; then
    if [[ "$machine" == URSA ]] ; then
	set +uxe
	echo "ERROR: Specify compiler when running rt.sh on Ursa." 1>&2
	echo "ERROR: Specify compiler: rt.sh [intel|intelllvm]" 1>&2
	exit 1
    else
	compiler=intel
    fi
fi

if [[ "$machine" == URSA ]] ; then
    runtime_log=$homedir/scripts/runtime.log.${machine}_${compiler}
else
    runtime_log=$homedir/scripts/runtime.log.$machine
fi

#set working directory
export workdir=${workdir:-"`pwd`/work-upp-${machine}-${compiler}"}
rm -rf $workdir
mkdir -p $workdir

#differentiates for orion and hercules
export rundir="${rundir}/upp-${machine}"
test -d "${rundir}" || mkdir -p "${rundir}"

#set log file
export rt_log=rt.log.${machine}_${compiler}
export logfile=`pwd`/$rt_log
if [ -f $logfile ] ; then
 rm -r $logfile
fi

#build executable
if [ "$build_exe" = "yes" ]; then
  cd ${test_v}
  mkdir -p ${test_v}/exec
  cd ${test_v}/tests
  ./compile_upp.sh -o upp_no_ifi.x -c "$compiler"
  status=$?
  if [ $status -eq 0 ]; then
    msg="Building executable successfully"
  else
    msg="Building executable with failure"
    postmsg "$logfile" "$msg"
    exit 2
  fi

  if [[ "$have_ifi" == yes && "$disable_ifi" == no ]] ; then
    ./compile_upp.sh -a -o upp_with_ifi.x -I -B -c "$compiler"
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

cat << EOF > $rt_log.temp
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
    echo "Warning: some tests exited with non-zero. status" >> $rt_log.temp
    echo >> $rt_log.temp
fi

cat "$rt_log" | grep "test:" >> $rt_log.temp
cat "$rt_log" | grep "baseline" >> $rt_log.temp
python ${test_v}/ci/rt-status.py >> $rt_log.temp
echo "===== End of UPP Regression Testing Log =====" >> $rt_log.temp
mv $rt_log.temp $rt_log
mv $rt_log ${test_v}/tests/logs

# should indicate failure to Jenkins
if [ $test_results -ne 0 ]; then
   python ${test_v}/ci/rt-status.py > changed_results.txt
   if [ $some_failed = YES ]; then
     echo "Warning: some tests exited with non-zero status." >> changed_results.txt
   fi
   exit 1
fi
