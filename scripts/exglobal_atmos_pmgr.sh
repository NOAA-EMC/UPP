#! /bin/ksh
#
#  Script name:         exglobal_atmos_pmgr.sh
#
#  This script monitors the progress of the gfs_post job
#  Script history log:
#  2021-05-25   Lin Gan: Ensure loganl.txt file created by atmos_analysis_calc before release atmos_post_anl job.
#
set -x

hour=00
typeset -Z2 hour

case $RUN in
   gfs) TEND=384
        TCP=385;;
  gdas) TEND=9
        TCP=10;; 
esac

if [ -e posthours ]; then
   rm -f posthours
fi

while [ $hour -lt $TCP ]; 
do
  echo $hour >>posthours
  if [ $hour -lt 120 ]
  then
     if [ $hour -eq 99 ]
     then
       typeset -Z3 hour
     fi
     let "hour=hour+1"
  else
     let "hour=hour+3"
  fi
done
postjobs=`cat posthours`

#
# Wait for all fcst hours to finish 
#
icnt=1
while [ $icnt -lt 1000 ]
do
  for fhr in $postjobs
  do 
    fhr3=`printf "%03d" $fhr`   
    if [ -s ${COMIN}/${RUN}.${cycle}.logf${fhr}.txt -o  -s ${COMIN}/${RUN}.${cycle}.logf${fhr3}.txt ]
    then
      if [ $fhr -eq 0 ]
      then 
        if [ -s ${COMIN}/${RUN}.${cycle}.loganl.txt ]
        then
          ecflow_client --event release_postanl
          ecflow_client --event release_post000
          postjobs=`echo $postjobs | sed "s/00//"`
        fi
      else    
        ecflow_client --event release_post${fhr3}
        # Remove current fhr from list
        postjobs=`echo $postjobs | sed "s/${fhr}//"`
      fi
    fi
  done
  
  result_check=`echo $postjobs | wc -w`
  if [ $result_check -eq 0 ]
  then
     break
  fi

  sleep 10
  icnt=$((icnt + 1))
  if [ $icnt -ge 1080 ]
  then
    msg="ABORTING after 3 hours of waiting for ${RUN} FCST hours $postjobs."
    err_exit $msg
  fi

done

echo Exiting $0

exit
