#!/bin/sh

#SBATCH -o out.fv3gefs_pe_test
#SBATCH -e out.fv3gefs_pe_test
#SBATCH -J fv3gefs_pe_test
#SBATCH -t 00:30:00
#SBATCH -N 4 --ntasks-per-node=12
##SBATCH -q debug
#SBATCH -q batch
#SBATCH -A ovp

set -x

# specify computation resource
export threads=1
export MP_LABELIO=yes
export OMP_NUM_THREADS=$threads
export APRUN="srun"

############################################
# Loading module
############################################
module purge
module use /contrib/spack-stack/spack-stack-1.8.0/envs/ue-intel-2021.5.0/install/modulefiles/Core
module load stack-intel/2021.5.0
module load stack-intel-oneapi-mpi/2021.5.1
module load libpng/1.6.37
module load jasper/2.0.32
module load prod_util/2.1.1
module load crtm/2.4.0.1
module list

msg="Starting fv3gefs test"
postmsg "$logfile" "$msg"

export cmp_grib2_grib2=/home/Wen.Meng/bin/cmp_grib2_grib2_new
export POSTGPEXEC=${svndir}/exec/upp.x

# specify forecast start time and hour for running your post job
export startdate=2022042400
export fhr=060
export CC=`echo $startdate | cut -c9-10`

# specify your running and output directory
export DATA=$rundir/fv3gefs_${startdate}_pe_test
rm -rf $DATA; mkdir -p $DATA
cd $DATA

export NEWDATE=`${NDATE} +${fhr} $startdate`
                                                                                       
export YY=`echo $NEWDATE | cut -c1-4`
export MM=`echo $NEWDATE | cut -c5-6`
export DD=`echo $NEWDATE | cut -c7-8`
export HH=`echo $NEWDATE | cut -c9-10`


cat > itag <<EOF
&model_inputs
fileName='$homedir/data_in/gefs/geaer.t${CC}z.atmf${fhr}.nemsio'
IOFORM='binarynemsiompiio'
grib='grib2'
DateStr='${YY}-${MM}-${DD}_${HH}:00:00'
MODELNAME='GFS'
fileNameFlux='$homedir/data_in/gefs/geaer.t${CC}z.sfcf${fhr}.nemsio'
/
 &NAMPGB
 KPO=47,PO=1000.,975.,950.,925.,900.,875.,850.,825.,800.,775.,750.,725.,700.,675.,650.,625.,600.,575.,550.,525.,500.,475.,450.,425.,400.,375.,350.,325.,300.,275.,250.,225.,200.,175.,150.,125.,100.,70.,50.,30.,20.,10.,7.,5.,3.,2.,1.,gocart_on=.true.,
/
EOF


rm -f fort.*

cp ${svndir}/fix/nam_micro_lookup.dat ./eta_micro_lookup.dat
cp $homedir/fix/postxconfig-NT-GEFS-CHEM.txt ./postxconfig-NT.txt

# copy flat files instead
#ens_pert_type=pos_pert_fcst
#sed < ${svndir}/parm/postxconfig-NT-GEFS.txt -e "s#negatively_pert_fcst#${ens_pert_type}#" > ./postxconfig-NT.txt

cp ${svndir}/parm/params_grib2_tbl_new ./params_grib2_tbl_new

cp ${svndir}/fix/chem/optics_luts_DUST.dat ./optics_luts_DUST.dat
cp ${svndir}/fix/chem/optics_luts_SALT.dat ./optics_luts_SALT.dat
cp ${svndir}/fix/chem/optics_luts_SOOT.dat ./optics_luts_SOOT.dat
cp ${svndir}/fix/chem/optics_luts_SUSO.dat ./optics_luts_SUSO.dat
cp ${svndir}/fix/chem/optics_luts_WASO.dat ./optics_luts_WASO.dat

export PGBOUT=pgbfile
${APRUN} ${POSTGPEXEC} < itag > outpost_nems_${NEWDATE}

#$COPYGB2 -x -i'4,0,80' -k'1 3 0 7*-9999 101 0 0' $PGBOUT tfile
#$WGRIB2 tfile -set_byte 4 11 1 -grib prmsl
#$COPYGB2 -x -i'4,1,5' -k'1 3 5 7*-9999 100 0 50000' $PGBOUT tfile
#$WGRIB2 tfile -set_byte 4 11 193 -grib h5wav
#cat  prmsl h5wav >> $PGBOUT
mv $PGBOUT geaer.t${CC}z.master.grb2f${fhr}

fhr2=`printf "%02d" $fhr`

filelist="geaer.t${CC}z.master.grb2f${fhr}"

for file in $filelist; do
export filein2=$file
ls -l ${filein2}
export err=$?

if [ $err = "0" ] ; then

 # use cmp to see if new pgb files are identical to the control one
 cmp ${filein2} $homedir/data_out/gefs/${filein2}.${machine}

 # if not bit-identical, use cmp_grib2_grib2 to compare each grib record
 export err1=$?
 if [ $err1 -eq 0 ] ; then
  msg="fv3gefs pe test: your new post executable generates bit-identical ${filein2} as the trunk"
  echo $msg
 else
  msg="fv3gefs pe test: your new post executable did not generate bit-identical ${filein2} as the trunk"
  echo $msg
  echo " start comparing each grib record and write the comparison result to *diff files"
  echo " check these *diff files to make sure your new post only change variables which you intend to change"
  $cmp_grib2_grib2 $homedir/data_out/gefs/${filein2}.${machine} ${filein2} > ${filein2}.diff
 fi

else

 msg="fv3gefs pe test: post failed using your new post executable to generate ${filein2}"
 echo $msg 2>&1 | tee -a TEST_ERROR

fi
postmsg "$logfile" "$msg"
done

echo "PROGRAM IS COMPLETE!!!!!" 2>&1 | tee SUCCESS
msg="Ending fv3gefs pe test"
postmsg "$logfile" "$msg"
