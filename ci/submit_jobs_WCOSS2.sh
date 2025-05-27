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
job_id=`qsub run_post_${test}_${machine}.sh`
export jobid_list=$jobid_list" "${job_id}
done

