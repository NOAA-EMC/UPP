#!/bin/bash
##########################################################################
# This script is used to submit test jobs on R&D machines.
# # Wen Meng, 05/2025, First version.
# ##########################################################################

export jobid_list=""

cd $workdir
for test in ${test_list}
do
  cp $svndir/ci/jobs-dev/run_post_${test}_${machine}.sh .
  job_id=$(sbatch --parsable -A "${accnr}" run_post_${test}_${machine}.sh)
  jobid_list="${jobid_list} ${job_id}"
done

#Run additional ifi tests
if [[ "$have_ifi" == "yes" && "$disable_ifi" == "no" ]] ; then
  for model in hrrr rrfs; do
    cp "$svndir/ci/jobs-dev/run_post_${model}_ifi_${machine}.sh" .
    job_id=$(sbatch --parsable -A "$accnr" "run_post_${model}_ifi_${machine}.sh")
    jobid_list="${jobid_list} $job_id" 
  
    cp "$svndir/ci/jobs-dev/run_ifi_standalone_${model}_${machine}.sh" .
    dep_job_id="$job_id"
    job_id=$(sbatch --parsable -A "$accnr" --dependency=afterany:"$dep_job_id" \
             "run_ifi_standalone_${model}_${machine}.sh")
    jobid_list="${jobid_list} $job_id"

    test_list="${test_list} ${model}_ifi"
  done

fi

export jobid_list
export test_list
