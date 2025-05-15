#!/bin/bash
##########################################################################
# This script is used to retrive runtime
# Wen Meng, 05/2025, First version.
##########################################################################

#get runtime for each test
export some_failed=NO
sleep 30
for job_id in $jobid_list; do
  ic=1
  sleep_loop_max=300
  while [ $ic -le $sleep_loop_max ]; do
     status=`qstat -x ${job_id} | awk 'FNR == 3' |  awk '{print $5}'`
     if [ "$status" = "F" ]; then
       break
     else
      ic=`expr $ic + 1`
      sleep 15
     fi
  done
  if [ $ic -lt $sleep_loop_max ]; then
     stime=`qstat -xf ${job_id} | grep stime | awk -F "=" '{print $2}'`
     stime=`date -d"$stime" +%s`
     etime=`qstat -xf ${job_id} | grep mtime | awk -F "=" '{print $2}'`
     etime=`date --date="$etime" +%s`
     runtime=$(( ($etime - $stime) ))
     runtime=`date -d@$runtime +%H:%M:%S`
     #runtime=`qstat -x ${job_id} | awk 'FNR == 3' | awk '{print $4}'`
     jobname=`qstat -x ${job_id} | awk 'FNR == 3' | awk '{print $2}'`
     runtime_b=`grep ${jobname} ${runtime_log} | awk '{print $2}' `
     echo "$runtime   $jobname ${runtime_b}"
     msg="Runtime: $jobname $runtime -- baseline ${runtime_b}"
     postmsg "$logfile" "$msg"
  fi
done

elapsed_time=$( printf '%02dh:%02dm:%02ds\n' $((SECONDS%86400/3600)) $((SECONDS%3600/60)) $((SECONDS%60)) )

python ${test_v}/ci/rt-status_${machine}.py
test_results=$?

if [ $some_failed = YES ] ; then
	  test_results=99
	    echo WARNING: some tests exited with non-zero status.
fi
