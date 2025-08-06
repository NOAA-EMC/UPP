#!/bin/bash

# Script to search for unclaimed index numbers (<post_avblfldidx>) in post_avblflds.xml
max_id=1015

i=1
while [ $i -le ${max_id} ]; do
   num_found=0
   while read -r line; do
      idx=`echo ${line} | grep -o -E [0-9]+`
      if [ ${idx} == $i ]; then
         num_found=1
         break
      fi
   done < out
   if [ ${num_found} -eq 0 ]; then
      echo "index number ${i} not in list"
   fi
   ((i++))
done
