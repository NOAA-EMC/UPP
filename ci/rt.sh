#!/bin/bash
######################################################################
# This script is desined for UPP regression tests run by UPP developer.
# Wen Meng, 12/2020, First version.
# Fernando Andrade-Maldonado 5/2023 rework for CLI Options
# Fernando Andrade-Maldonado / Wen Meng 9/2023 Add Hercules, fix typos, and refactor
# Fernando Andrade-Maldonado 4/2024 Additional Log info
# Wen Meng 05/2025 Refactor to support WCOSS2 and R&D machines
######################################################################
set -xue
SECONDS=0

git_branch="develop"
git_url="https://github.com/NOAA-EMC/UPP.git"
clone_on="no"
export disable_ifi="no" # don't use libIFI, even if it is present
build_exe="yes" #build executable

while getopts a:w:h:r:t:b:u:cde opt; do
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
    e) build_exe="no" # don't build executable
        ;;
  esac
done

#UPP working copy
export test_v=${test_v:-`pwd`/..}
if [[ $clone_on == "yes" ]]; then
  rm -rf $test_v
  mkdir -p $test_v
  git clone -b $git_branch $git_url $test_v
fi
export svndir=${test_v}

if [[ -d $svndir/sorc/libIFI.fd/src/ ]] ; then
    export have_ifi=yes
else
    export have_ifi=no
fi

#find machine
mac=$(hostname | cut -c1-1)
mac2=$(hostname | cut -c1-2)
mac3=$(hostname | cut -c1-4)
if [ $mac2 = hf ]; then # for HERA
 export machine=HERA
 export homedir=${homedir:-"/scratch2/NAGAPE/epic/UPP/test_suite"}
 export rundir=${rundir:-"/scratch1/NCEPDEV/stmp2/${USER}"}
 export accnr=${accnr:-"rtrr"}
 module purge
 module use /contrib/spack-stack/spack-stack-1.8.0/envs/ue-intel-2021.5.0/install/modulefiles/Core
 module load stack-intel/2021.5.0
 module load stack-intel-oneapi-mpi/2021.5.1
 module load prod_util/2.1.1
elif [ $mac3 = orio ] ; then
 export machine=ORION
 export homedir=${homedir:-"/work/noaa/epic/UPP"}
 export rundir=${rundir:-"/work2/noaa/stmp/$USER"}
 export accnr=${accnr:-"rtrr"}
 module purge
 module use /apps/contrib/spack-stack/spack-stack-1.8.0/envs/ue-intel-2021.9.0/install/modulefiles/Core
 module load stack-intel/2021.9.0
 module load stack-intel-oneapi-mpi/2021.9.0
 module load prod_util/2.1.1
 module load python/3.10.8
elif [ $mac3 = herc ] ; then
 export machine=HERCULES
 export homedir=${homedir:-"/work/noaa/epic/UPP"}
 export rundir=${rundir:-"/work2/noaa/stmp/$USER"}
 export accnr=${accnr:-"rtrr"}
 module purge
 module use /apps/contrib/spack-stack/spack-stack-1.8.0/envs/ue-intel-2021.9.0/install/modulefiles/Core
 module load stack-intel/2021.9.0
 module load stack-intel-oneapi-mpi/2021.9.0
 module load prod_util/2.1.1
 module load python/3.10.8
elif [ $mac = d -o $mac = c ]; then #for WCOSS2
 export machine=WCOSS2
 export homedir=${homedir:-"/u/wen.meng/noscrub/ncep_post/post_regression_test_new"}
 export rundir=${rundir:-"/lfs/h2/emc/ptmp/$USER"}
 export accnr=${accnr:-"GFS-DEV"}
 module reset
 module load intel/19.1.3.304
 module load PrgEnv-intel/8.1.0
 module load craype/2.7.8
 module load cray-mpich/8.1.7
 module load prod_util/2.0.14
 module load python/3.12.0
fi

#set working directory
export workdir=${workdir:-"`pwd`/work-upp-${machine}"}
rm -rf $workdir
mkdir -p $workdir

#differentiates for orion and hercules
export rundir="${rundir}/upp-${machine}"
#test -d "${rundir}" || mkdir -p "${rundir}"
rm -rf ${rundir}; mkdir -p ${rundir}

#set log file
export logfile=`pwd`/rt.log.$machine
if [ -f $logfile ] ; then
 rm -r $logfile
fi
export runtime_log=$svndir/ci/runtime.log.$machine

#build executable
if [ "$build_exe" == "yes" ]; then
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

  if [[ "$have_ifi" == "yes" && "$disable_ifi" == "no" ]] ; then
    if [[ "${machine}" == "WCOSS2" ]]; then ##No ifi standalone executable
      ./compile_upp.sh -a -o upp_with_ifi.x -I 
      status=$?
    else
      ./compile_upp.sh -a -o upp_with_ifi.x -I -B
      status=$?
    fi
    if [ "$status" -eq 0 ]; then
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

#Setting tests
export test_list="nmmb fv3gefs fv3r fv3r_ifi_missing hrrr rap fv3hafs 3drtma fv3gfs"

#submit test jobs
cd $svndir/ci
if [ "${machine}" = "WCOSS2" ]; then
  source "./submit_jobs_${machine}.sh"
else  ##R&D machines
  source "./submit_jobs.sh"
fi

set +xe
echo "Job cards submitted for enabled tests, waiting on timestamps for finished jobs..."

#get run time for each test
cd $svndir/ci
if [ "${machine}" = "WCOSS2" ]; then
  source "./check_runtime_${machine}.sh"
else
  source "./check_runtime.sh"
fi
