#!/bin/bash

run=${1:-gfs}
bdate=${2:-2018100400}
edate=${2:-2018100400}
export allfhr=${3:-"030"}

#Input Data
#export COMINP=/u/Wen.Meng/ptmp/prfv3rt1
export COMINP=/scratch3/NCEPDEV/ovp/Wen.Meng/prfv3test
#export COMINP=/scratch3/NCEPDEV/stmp1/Alexei.A.Belochitski/tmp/comrot/test4

#Working directory
tmp=/scratch3/NCEPDEV/stmp2/$USER/nceppost
mkdir -p $tmp/ecf
diroutp=$tmp/outputs
mkdir -p $diroutp

#UPP location
export svndir=`pwd`

while [[ $bdate -le $edate ]]; do
   yyyymmdd=`echo $bdate | cut -c1-8`
   sed -e "s|CURRENTDATE|$bdate|" \
       -e "s|STDDIR|$diroutp|" \
       -e "s|RRR|$run|" \
      run_JGLOBAL_NCEPPOST_theia-new >$tmp/ecf/run_JGLOBAL_NCEPPOST_${run}.$bdate
   sbatch $tmp/ecf/run_JGLOBAL_NCEPPOST_${run}.$bdate
   bdate=`ndate +24 $bdate`
done
