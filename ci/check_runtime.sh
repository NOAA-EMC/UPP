#!/bin/bash
##########################################################################
# This script is used to retrive runtime on R&D machines.
# Wen Meng, 05/2025, First version.
##########################################################################

#get runtime for each test
export some_failed="NO"
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

if [ "$some_failed" = "YES" ] ; then
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
   if [ "$some_failed" = "YES" ]; then
     echo "Warning: some tests exited with non-zero status." >> changed_results.txt
   fi
   exit 1
fi
