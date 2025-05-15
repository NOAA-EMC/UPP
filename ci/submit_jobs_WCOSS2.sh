#!/bin/bash

jobid_list=""
#test_list="hrrr"
test_list="nmmb_Grib2 fv3gefs fv3r hrrr rap fv3hafs 3drtma fv3gfs"

cd $workdir
for test in ${test_list}
do
cp $svndir/ci/jobs-dev/run_post_${test}_${machine}.sh .
job_id=`qsub run_post_${test}_${machine}.sh`
export jobid_list=$jobid_list" "${job_id}
done
