#!/bin/bash
##########################################################################
# This script is used to submit test jobs on WCOSS2.
# # Wen Meng, 05/2025, First version.
# ##########################################################################

jobid_list=""

cd $workdir
for test in ${test_list}
do
cp $svndir/ci/jobs-dev/run_post_${test}_${machine}.sh .
job_id=`qsub -A ${accnr} run_post_${test}_${machine}.sh`
export jobid_list=$jobid_list" "${job_id}
done

#Run additional ifi tests
if [[ "$have_ifi" == yes && "$disable_ifi" == no ]] ; then
  cp $svndir/ci/jobs-dev/run_post_hrrr_ifi_${machine}.sh .
  job_id=`qsub -A ${accnr} run_post_hrrr_ifi_${machine}.sh`
  export jobid_list=$jobid_list" "${job_id}
  export test_list=${test_list}" hrrr_ifi"

  cp $svndir/ci/jobs-dev/run_post_fv3r_ifi_${machine}.sh .
  job_id=`qsub -A ${accnr} run_post_fv3r_ifi_${machine}.sh`
  export jobid_list=$jobid_list" "${job_id}
  export test_list=${test_list}" fv3r_ifi"
fi
