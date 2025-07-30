#!/bin/sh 
 
#PBS -o out.3drtma.pe_test
#PBS -e out.3drtma.pe_test
#PBS -N 3drtma.pe_test
#PBS -l walltime=00:30:00
#PBS -q debug
#PBS -A RTMA-DEV
#PBS -l place=vscatter,select=1:ncpus=72
#PBS -V

set -x

# specify computation resource
export threads=1
export OMP_NUM_THREADS=$threads
export APRUN="mpiexec -l -n 72 -ppn 72"

echo "starting time"
date

############################################
# Loading module
############################################
module reset
module load intel/19.1.3.304
module load PrgEnv-intel/8.1.0
module load craype/2.7.8
module load cray-mpich/8.1.7
module load cray-pals/1.0.12
module load hdf5/1.10.6
module load netcdf/4.7.4
module load libjpeg/9c
module load prod_util/2.0.8
module list

msg="Starting rtma pe test"
postmsg "$logfile" "$msg"

export cmp_grib2_grib2=/u/wen.meng/bin/cmp_grib2_grib2_new
export POSTGPEXEC=${svndir}/exec/upp.x

# specify your running and output directory
export startdate=2025071400
export DATA=$rundir/rtma_${startdate}_pe_test

export NEWDATE=$startdate

export YY=`echo ${NEWDATE} | cut -c1-4`
export MM=`echo ${NEWDATE} | cut -c5-6`
export DD=`echo ${NEWDATE} | cut -c7-8`
export HH=`echo ${NEWDATE} | cut -c9-10`
export min=00

rm -rf $DATA; mkdir -p $DATA
cd $DATA

cat > itag <<EOF
&model_inputs
fileName='$homedir/data_in/3drtma/rtma3d.t${HH}z.wrf_inout.nc'
IOFORM='netcdf'
grib='grib2'
DateStr='${YY}-${MM}-${DD}_${HH}:${min}:00'
MODELNAME='RAPR'
SUBMODELNAME='RTMA'
/
&NAMPGB
KPO=47,PO=2.,5.,7.,10.,20.,30.,50.,70.,75.,100.,125.,150.,175.,200.,225.,250.,275.,300.,325.,350.,375.,400.,425.,450.,475.,500.,525.,550.,575.,600.,625.,650.,675.,700.,725.,750.,775.,800.,825.,850.,875.,900.,925.,950.,975.,1000.,1013.2
/
EOF

#copy fix data
cp $homedir/fix/fix_2.3.0/*bin .
cp ${svndir}/parm/params_grib2_tbl_new params_grib2_tbl_new
cp ${svndir}/parm/postxconfig-NT-3drtma.txt postxconfig-NT.txt
cp ${svndir}/fix/rap_micro_lookup.dat eta_micro_lookup.dat

${APRUN} ${POSTGPEXEC} < itag > wrfpost2.out

# operational rtma post processing generates 3 files
filelist="WRFTWO.GrbF00 \
          WRFPRS.GrbF00 \
          WRFNAT.GrbF00"

for file in $filelist; do
export filein2=$file
ls -l ${filein2}
export err=$?

if [ $err = "0" ] ; then

 # operational rtma post processing generates 3 files, start with BGDAWP first
 # use cmp to see if new pgb files are identical to the control one
 cmp ${filein2} $homedir/data_out/${filein2}.${machine}

 # if not bit-identical, use cmp_grib2_grib2 to compare each grib record
 export err1=$?
 if [ $err1 -eq 0 ] ; then
  msg="rtma pe test: your new post executable generates bit-identical ${filein2} as the trunk"
  echo $msg
 else
  #msg="rtma pe test: your new post executable did not generate bit-identical ${filein2} as the trunk"
  #echo $msg
  echo " start comparing each grib record and write the comparison result to *diff files"
  echo " check these *diff files to make sure your new post only change variables which you intend to change"
  $cmp_grib2_grib2 $homedir/data_out/${filein2}.${machine} ${filein2} > ${filein2}.diff
  cmp ${filein2}.diff $homedir/data_out/${filein2}.diff
  export err2=$?
  if [ $err2 -eq 0 ] ; then
   msg="rtma pe test: your new post executable is fine in ${filein2}"
   echo $msg
  else
   msg="rtma pe test: your new post executable did generate changed results in ${filein2}"
   echo $msg
  fi
 fi


else

 msg="rtma pe test: post failed using your new post executable to generate ${filein2}"
 echo $msg

fi
postmsg "$logfile" "$msg"
done

echo "PROGRAM IS COMPLETE!!!!!"
msg="Ending rap pe test"
postmsg "$logfile" "$msg"


