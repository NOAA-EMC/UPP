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
job_id=$(qsub -A "${accnr}" run_post_${test}_${machine}.sh)
jobid_list="${jobid_list} ${job_id}"
done

#Run additional ifi tests
if [[ "$have_ifi" == "yes" && "$disable_ifi" == "no" ]] ; then
  for ifi_test in hrrr_ifi fv3r_ifi; do
    cp $svndir/ci/jobs-dev/run_post_${ifi_test}_${machine}.sh .
    job_id=$(qsub -A "${accnr}" run_post_${ifi_test}_${machine}.sh)
    jobid_list="${jobid_list} ${job_id}"
    test_list=${test_list}" ${ifi_test}"
  done
fi

export jobid_list
export test_list
