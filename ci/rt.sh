#!/bin/bash
######################################################################
# This script is desined for UPP regression tests run by UPP developer.
# Wen Meng, 12/2020, First version.
# Fernando Andrade-Maldonado 5/2023 rework for CLI Options
# Fernando Andrade-Maldonado / Wen Meng 9/2023 Add Hercules, fix typos, and refactor
# Fernando Andrade-Maldonado 4/2024 Additional Log info
# Wen Meng 05/2025 Refactor to support WCOSS2 and R&D machines
# Sam Trahan 06/2025 Add usage message, Ursa support, and multi-compiler support
# Gillian Petro 06/2025 Update to spack-stack 1.9.1; require compiler indication on Orion/Hercules
# Wen Meng and Ben Blake, 07/2025, Update test names, add RRFS, MPAS, DAFS, SFS tests
######################################################################
set -xue
SECONDS=0

git_branch="develop"
git_url="https://github.com/NOAA-EMC/UPP.git"
clone_on="no"
disable_ifi="no" # don't use libIFI, even if it is present
disable_gtg="no" # don't use post_gtg, even if it is present
print_full_help="no"
build_exe="yes" #build executable
compiler="MISSING"

usage() {
  set +xue

  cat<<EOF

Usage: rt.sh -a account -C compiler -r /path/to/scrub/space [-options] [compiler]
Executes UPP regression tests. Includes IFI tests if ../sorc/libIFI.fd exists.
Includes GTG tests if ../sorc/ncep_post.fd/post_gtg.fd exists.

Results are here:
  ../tests/logs/MACHINE_compiler.log = report of regression tests for each machine and compiler.
  changed_results.txt = A list of tests whose results have changed.

Always set these:
  -a account = accounting code for job submission. Default account is often overused. Always set this!
  -C = chosen compiler. (Capital C) Default: intel. Mandatory on Ursa!
  -r rundir = path to a scrub space. Default area is often over quota. Always set this!

General options:
  -d = disable ifi tests even if ifi is available
  -g = disable gtg tests even if gtg is available
  -e = don't build the UPP executable
  -h homedir = path to the regression test data
  -w workdir = directory to store per-job batch and log files.
  -l test_list = list of RTs to run in a string (e.g., "sfs gefsv12 rap" )
EOF

  if [[ "$print_full_help" == YES ]] ; then
    cat<<EOF
  -H = print this message and exit.

Special run mode: run rt.sh outside the repository. Automatically clones the repository.
Syntax: rt.sh -a account -r /path/to/scrub/space -c -u url -b branch [options] [compiler]

Additional options:
  -c = Tells rt.sh it is running outside a repository. (Lower-case c)
  -t test_v = Location to clone the repository. Default: Overwrite .. with the clone.
  -u url = Mandatory: URL of a repository to clone. Not for general use.
  -b branch = Mandatory: branch in the repository to clone
EOF
  else
      cat<<EOF
  -H = print full help message and exit. Includes special-use options.
EOF
  fi

  if [[ "$#" -gt 0 ]] ; then
    echo
    echo "------------------------------------------------------------------------"
    echo "$@"
    echo "------------------------------------------------------------------------"
  fi
}

check_for_dash() {
  if [[ -z "${OPTARG}" ]] ; then
    echo "Argument error: -$opt argument is the empty string"
    usage FATAL ERROR: Script is exiting due to invalid argument. See error message above. 1>&2
    exit 2
  fi
  if [[ "${OPTARG:0:1}" == '-' ]] ; then
    echo "Argument error: -$opt requires an argument"
    usage FATAL ERROR: Script is exiting due to a missing argument. See error message above 1>&2
    exit 2
  fi
}

# Space required at start and and of string for pattern matching in check_valid_tests
valid_tests=' sfs gefsv12 gefsv13 nmmb rap hrrr hafs 3drtma mpas mpas_hfip rrfs rrfs_ifi_missing gfs '

check_valid_tests() {
   local tests=${@}
   if [[ -n ${tests} ]]; then
      test_list=''
      read -a tests_to_run <<< ${tests}
      for t in ${tests_to_run[@]}
      do
         if [[ ${valid_tests} =~ " ${t} " ]]; then
            test_list+="${t} "
         else
            echo "${t} is not a valid test"
         fi
      done
      if [[ -z ${test_list} ]]; then
         echo "No valid tests provided. Exiting..."
	 exit 1
      fi
      export test_list=${test_list}
   fi
   echo "rt.sh will run ${test_list}"
}

set +x
export OPTERR=1
while getopts a:w:h:r:l:t:b:u:C:cdgHe opt; do
  case $opt in
    C) compiler=${OPTARG} ; check_for_dash
        ;;
    d) disable_ifi=yes
        ;;
    g) disable_gtg=yes
        ;;
    a) accnr=${OPTARG} ; check_for_dash
        ;;
    w) workdir=${OPTARG} ; check_for_dash
        ;;
    h) homedir=${OPTARG} ; check_for_dash
        ;;
    r) rundir=${OPTARG} ; check_for_dash
        ;;
    l) check_valid_tests ${OPTARG}
       ;;
    t) test_v=${OPTARG} ; check_for_dash
        ;;
    b) git_branch=${OPTARG} ; check_for_dash
        ;;
    u) git_url=${OPTARG} ; check_for_dash
        ;;
    c) clone_on="yes"
        ;;
    e) build_exe="no" # don't build executable
        ;;
    H) print_full_help=YES ; usage ; exit 1
        ;;
    :) echo "Error: Required argument not provided." 
       help 
       exit 2 
       ;;
    *)
       usage FATAL ERROR: Invalid -option. See error message above. 1>&2
       exit 2
        ;;
  esac
done

# Fail if positional arguments are present:
positional_count=$(( $# - OPTIND + 1 ))
if (( positional_count > 0)) ; then
  if (( positional_count > 1)) ; then
    arguments=arguments
  else
    arguments=argument
  fi
  shift $(( OPTIND - 1 ))
  usage FATAL ERROR: Script is aborting due to spurious $arguments: "$@" 2>&1
  exit 2
fi
set -x

# Set test list if not set

test_list=${test_list:-${valid_tests}}

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

if [[ -f $svndir/sorc/ncep_post.fd/post_gtg.fd/gtg.config.hrrr ]] ; then
    export have_gtg=yes
else
    export have_gtg=no
fi

#find machine
mac=$(hostname | cut -c1-1)
mac2=$(hostname | cut -c1-2)
mac3=$(hostname | cut -c1-4)
if [ $mac2 = hf ]; then # for HERA
 export machine=HERA
 export homedir=${homedir:-"/scratch4/NAGAPE/epic/role-epic/hera/UPP_test_suite"}
 export rundir=${rundir:-"/scratch3/NCEPDEV/stmp/${USER}"}
 export accnr=${accnr:-"rtrr"}
 module purge
 module use /contrib/spack-stack/spack-stack-1.9.2/envs/ue-oneapi-2024.2.1/install/modulefiles/Core
 module use /contrib/spack-stack/spack-stack-1.9.2/envs/ue-oneapi-2024.2.1/install/modulefiles/intel-oneapi-mpi/2021.13-sbi3u54/gcc/13.3.0
 module load stack-oneapi/2024.2.1
 module load stack-intel-oneapi-mpi/2021.13
 module load prod_util/2.1.1
elif [ $mac2 = uf ]; then # for Ursa
 export machine=URSA
 export homedir=${homedir:-"/scratch4/NAGAPE/epic/role-epic/ursa/UPP/test_suite"}
 export rundir=${rundir:-"/scratch3/NCEPDEV/stmp/$USER/scrub"}
 export accnr=${accnr:-"rtrr"}
 module use /contrib/spack-stack/spack-stack-1.9.2/envs/ue-oneapi-2024.2.1/install/modulefiles/Core
 module use /contrib/spack-stack/spack-stack-1.9.2/envs/ue-oneapi-2024.2.1/install/modulefiles/intel-oneapi-mpi/2021.13-haww6b3/gcc/12.4.0
 module load stack-oneapi/2024.2.1
 module load stack-intel-oneapi-mpi/2021.13
 module load prod_util/2.1.1
 module load python/3.11.7
elif [ $mac3 = orio ] ; then
 export machine=ORION
 export homedir=${homedir:-"/work/noaa/epic/role-epic/orion/UPP"}
 export rundir=${rundir:-"/work2/noaa/stmp/$USER"}
 export accnr=${accnr:-"rtrr"}
 module purge
 module use /apps/contrib/spack-stack/spack-stack-1.9.2/envs/ue-oneapi-2024.1.0/install/modulefiles/Core
 module use /apps/contrib/spack-stack/spack-stack-1.9.2/envs/ue-oneapi-2024.1.0/install/modulefiles/intel-oneapi-mpi/2021.13-li242lf/gcc/12.2.0
 module load stack-oneapi/2024.2.1
 module load stack-intel-oneapi-mpi/2021.13
 module load prod_util/2.1.1
 module load python/3.11.7
elif [ $mac3 = herc ] ; then
 export machine=HERCULES
 export homedir=${homedir:-"/work/noaa/epic/role-epic/hercules/UPP"}
 export rundir=${rundir:-"/work2/noaa/stmp/$USER"}
 export accnr=${accnr:-"rtrr"}
 module purge
 module use /apps/contrib/spack-stack/spack-stack-1.9.2/envs/ue-oneapi-2024.1.0/install/modulefiles/Core
 module use /apps/contrib/spack-stack/spack-stack-1.9.2/envs/ue-oneapi-2024.1.0/install/modulefiles/intel-oneapi-mpi/2021.13-sqiixt7/gcc/13.3.0
 module load stack-oneapi/2024.2.1
 module load stack-intel-oneapi-mpi/2021.13
 module load prod_util/2.1.1
 module load python/3.11.7
elif [ $mac = d -o $mac = c ]; then #for WCOSS2
 export machine=WCOSS2
 export homedir=${homedir:-"/lfs/h2/emc/vpppg/noscrub/wen.meng/test_suite"}
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

if [[ "$compiler" == MISSING ]] ; then
   if [[ "$machine" == "URSA" ]]; then
	   usage FATAL ERROR: You must specify the compiler on Ursa: -C 'intel|intelllvm' 1>&2
	   exit 2
   else
	   compiler=intel
   fi
fi

export compiler

#set working directory
export workdir=${workdir:-"`pwd`/work-upp-${machine}-${compiler}"}
rm -rf $workdir
mkdir -p $workdir

export cmp_grib2_grib2=$svndir/ci/cmp_grib2_grib2.sh

#differentiates for orion and hercules
export rundir="${rundir}/upp-${machine}"
#test -d "${rundir}" || mkdir -p "${rundir}"
rm -rf ${rundir}; mkdir -p ${rundir}

#set log file
export rt_log=rt.log.${machine}_${compiler}
export logfile=`pwd`/$rt_log
if [ -f $logfile ] ; then
 rm -r $logfile
fi
export runtime_log=$svndir/ci/runtime.log.${machine}_${compiler}

# Validate post_avblflds.xml against schema
cd ${test_v}/parm
if [[ ${machine} != "URSA" ]]; then
  xmllint --noout --schema EMC_POST_Avblflds_Schema.xsd post_avblflds.xml
  # If an error results from running xmllint, rt.sh will terminate execution and report that the XML fails to validate. 
fi

#build executable
if [ "$build_exe" == "yes" ]; then
  cd ${test_v}
  mkdir -p ${test_v}/exec
  cd ${test_v}/tests
  ./compile_upp.sh -o upp_no_ifi_gtg.x -c "$compiler"
  status=$?
  if [ $status -eq 0 ]; then
    msg="Building executable successfully"
  else
    msg="Building executable with failure"
    postmsg "$logfile" "$msg"
    exit 2
  fi

  if [[ "$have_ifi" == "yes" && "$disable_ifi" == "no" && "$have_gtg" == "yes" && "$disable_gtg" == "no" ]] ; then
    if [[ "${machine}" == "WCOSS2" ]]; then ##GTG tests only supported on WCOSS2
      ./compile_upp.sh -a -o upp_with_ifi_gtg.x -I -g -c "$compiler"
      status=$?
    else
      msg="GTG tests are not currently supported on machines other than WCOSS2, exiting"
      postmsg "$logfile" "$msg"
      exit 2
    fi
    if [ "$status" -eq 0 ]; then
      msg="Built UPP+IFI+GTG executables successfully"
    else
      msg="Built UPP+IFI+GTG executables with failure"
      postmsg "$logfile" "$msg"
      exit 2
    fi
    ln -s upp_with_ifi_gtg.x $svndir/exec/upp.x
  elif [[ "$have_ifi" == "yes" && "$disable_ifi" == "no" && "$have_gtg" == "no" ]] ; then
    if [[ "${machine}" == "WCOSS2" ]]; then ##No ifi standalone executable
      ./compile_upp.sh -a -o upp_with_ifi.x -I -c "$compiler"
      status=$?
    else
      ./compile_upp.sh -a -o upp_with_ifi.x -I -B -c "$compiler"
      status=$?
    fi
    if [ "$status" -eq 0 ]; then
      msg="Built UPP+IFI executables successfully"
    else
      msg="Built UPP+IFI executables with failure"
      postmsg "$logfile" "$msg"
      exit 2
    fi
    ln -s upp_with_ifi.x $svndir/exec/upp.x
  elif [[ "$have_gtg" == "yes" && "$disable_gtg" == "no" && "$have_ifi" == "no" ]] ; then
    if [[ "${machine}" == "WCOSS2" ]]; then ##GTG tests only supported on WCOSS2
      ./compile_upp.sh -a -o upp_with_gtg.x -g -c "$compiler"
      status=$?
    else
      msg="GTG tests are not currently supported on machines other than WCOSS2, exiting"
      postmsg "$logfile" "$msg"
      exit 2
    fi
    if [ "$status" -eq 0 ]; then
      msg="Built UPP+GTG executables successfully"
    else
      msg="Built UPP+GTG executables with failure"
      postmsg "$logfile" "$msg"
      exit 2
    fi
    ln -s upp_with_gtg.x $svndir/exec/upp.x
  else
    ln -s upp_no_ifi_gtg.x $svndir/exec/upp.x
  fi

  postmsg "$logfile" "$msg"
fi

# Create job cards from template for RDHPCS
if [[ ${machine} != "WCOSS2" ]]; then
   
   cd $svndir/ci/jobs-dev

   source machine.sh
   source test.sh
   source atparse.sh

   for test in ${test_list}
   do
      set_global
      ${test}
      atparse < run_post_${test}_template.sh > run_post_${test}_${machine}.sh
   done

fi

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
