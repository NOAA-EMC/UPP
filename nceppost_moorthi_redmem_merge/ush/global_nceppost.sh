#!/bin/ksh
################################################################################
####  UNIX Script Documentation Block
#                      .                                             .
# Script name:         global_nceppost.sh           
# Script description:  Posts the global pressure GRIB file
#
# Author:        Mark Iredell       Org: NP23         Date: 1999-05-01
#
# Abstract: This script reads a single global GFS IO file and (optionally)
#   a global flux file and creates a global pressure GRIB file.
#   The resolution and generating code of the output GRIB file can also
#   be set in the argument list.
#
# Script history log:
# 1999-05-01  Mark Iredell
# 2007-04-04  Huiya Chuang: Modify the script to run unified post
# 2012-03-00  Shrinivas Moorthi - update for gaea and zeus
# 2012-06-04  Jun Wang: add grib2 option
# 2012-xx-xx  Shrinivas Moorthi - added multi-machine option and generalized defaults
# 2015-03-20  Lin Gan: add Perl for Post XML performance upgrade
# 2016-02-08  Lin Gan: Modify to use Vertical Structure
#
# Usage:  global_postgp.sh SIGINP FLXINP FLXIOUT PGBOUT PGIOUT IGEN
#
#   Input script positional parameters:
#     1             Input sigma file
#                   defaults to $SIGINP
#     2             Input flux file
#                   defaults to $FLXINP
#     3             Output flux index file
#                   defaults to $FLXIOUT
#     4             Output pressure GRIB file
#                   defaults to $PGBOUT
#     5             Output pressure GRIB index file
#                   defaults to $PGIOUT, then to none
#     8             Model generating code,
#                   defaults to $IGEN, then to input sigma generating code
#
#   Imported Shell Variables:
#     SIGINP        Input sigma file
#                   overridden by $1
#     FLXINP        Input flux file
#                   overridden by $2
#     FLXIOUT       Output flux index file
#                   overridden by $3
#     PGBOUT        Output pressure GRIB file
#                   overridden by $4. If not defined,
#                   post will use the filename specified in
#                   the control file
#     PGIOUT        Output pressure GRIB index file
#                   overridden by $5; defaults to none
#     IGEN          Model generating code
#                   overridden by $8; defaults to input sigma generating code
##### Moorthi: Add new imported shell variable for running chgres
#     CHGRESSH      optional: the script to run chgres
#		    default to ${USHGLOBAL}/global_chgres.sh
#     SIGLEVEL      optional: the coordinate text file
#		    no default - uses coordinate from input sigma file
##### Chuang: Add new imported Shell Variable for ncep post
#     OUTTYP        Output file type for chgres
#                   1: if user has a sigma file and needs post to run chgres to convert to gfs io file
#                   2: if user already has a gfs io file
#                   3: if user uses post to read sigma file directly
#                   4: if user uses post to read nemsio file directly
#                   0: if user wishes to generate both gfsio and sigma files
#     VDATE         Verifying date 10 digits yyyymmddhh
#     GFSOUT        Optional, output file name from chgres which is input file name to nceppost 
#                   if model already runs gfs io, make sure GFSOUT is linked to the gfsio file
#     CTLFILE       Optional, Your version of control file if not using operational one
#     OVERPARMEXEC  Optional, the executable for changing Grib KPDS ID
#		    default to to ${EXECGLOBAL}/overparm_grib 
#     CHGRESTHREAD  Optional, speed up chgres by using multiple threads 
#		    default to 1
#     FILTER        Optional, set to 1 to filter SLP and 500 mb height using copygb
#     D3DINP        Optional, Inout D3D file, if not defined, post will run
#                   without processing D3D file
#     D3DOUT        Optional, output D3D file, if not defined, post will
#                   use the file name specified in the control file
#     IPVOUT        Optional, output IPV file, if not defined, post will
#                   use the file name specified in the control file
#    GENPSICHI      Optional, set to YES will generate psi and chi and
#                   append it to the end of PGBOUT.  Default to NO
#    GENPSICHIEXE   Optional, specify where executable is for generating
#                   psi and chi.
########################################################################       
#     EXECUTIL      Directory for utility executables
#                   defaults to /nwprod/util/exec
#     USHUTIL       Directory for utility scripts
#                   defaults to /nwprod/util/ush
#     EXECGLOBAL    Directory for global executables
#                   defaults to /nwprod/exec
#     USHGLOBAL     Directory for global scripts
#                   defaults to /nwprod/ush
#     DATA          working directory
#                   (if nonexistent will be made, used and deleted)
#                   defaults to current working directory
#     MP            Multi-processing type ("p" or "s")
#                   defaults to "p", or "s" if LOADL_STEP_TYPE is not PARALLEL
#     XC            Suffix to add to executables
#                   defaults to none
#     POSTGPEXEC    Global post executable
#                   defaults to ${EXECGLOBAL}/ncep_post
#     GRBINDEX      GRIB index maker
#                   defaults to ${EXECUTIL}/grbindex$XC
#     ANOMCATSH     Global anomaly GRIB script
#                   defaults to ${USHGLOBAL}/global_anomcat.sh
#     POSTGPLIST    File containing further namelist inputs
#                   defaults to /dev/null
#     INISCRIPT     Preprocessing script
#                   defaults to none
#     LOGSCRIPT     Log posting script
#                   defaults to none
#     ERRSCRIPT     Error processing script
#                   defaults to 'eval [[ $err = 0 ]]'
#     ENDSCRIPT     Postprocessing script
#                   defaults to none
#     POSTGPVARS    Other namelist inputs to the global post executable
#                   such as IDRT,KO,PO,KTT,KT,PT,KZZ,ZZ,
#                           NCPUS,MXBIT,IDS,POB,POT,MOO,MOOA,MOW,MOWA,
#                           ICEN,ICEN2,IENST,IENSI
#                   defaults to none set
#     NTHREADS      Number of threads
#                   defaults to 1
#     NTHSTACK      Size of stack per thread
#                   defaults to 64000000
#     VERBOSE       Verbose flag (YES or NO)
#                   defaults to NO
#     PGMOUT        Executable standard output
#                   defaults to $pgmout, then to '&1'
#     PGMERR        Executable standard error
#                   defaults to $pgmerr, then to '&1'
#     pgmout        Executable standard output default
#     pgmerr        Executable standard error default
#     REDOUT        standard output redirect ('1>' or '1>>')
#                   defaults to '1>', or to '1>>' to append if $PGMOUT is a file
#     REDERR        standard error redirect ('2>' or '2>>')
#                   defaults to '2>', or to '2>>' to append if $PGMERR is a file
#
#   Exported Shell Variables:
#     PGM           Current program name
#     pgm
#     ERR           Last return code
#     err
#
#   Modules and files referenced:
#     scripts    : $INISCRIPT
#                  $LOGSCRIPT
#                  $ERRSCRIPT
#                  $ENDSCRIPT
#                  $ANOMCATSH
#
#     programs   : $POSTGPEXEC
#                  $GRBINDEX
#
#     input data : $1 or $SIGINP
#                  $2 or $SFCINP
#                  $POSTGPLIST
#
#     output data: $3 or $FLXIOUT
#                  $4 or $PGBOUT
#                  $5 or $PGIOUT
#                  $PGMOUT
#                  $PGMERR
#
#     scratch    : ${DATA}/postgp.inp.sig
#                  ${DATA}/postgp.inp.flx
#                  ${DATA}/postgp.out.pgb
#
# Remarks:
#
#   Condition codes
#      0 - no problem encountered
#     >0 - some problem encountered
#
#   Control variable resolution priority
#     1 Command line argument.
#     2 Environment variable.
#     3 Inline default.
#
# Attributes:
#   Language: POSIX shell
#   Machine: WCOSS, ZEUS, THEIA, GAEA
#
####
################################################################################
#  Set environment.
export VERBOSE=${VERBOSE:-"NO"}
if [[ "$VERBOSE" = "YES" ]] ; then
   echo $(date) EXECUTING $0 $* >&2
   set -x
fi
export machine=${machine:-WCOSS}
export machine=$(echo $machine|tr '[a-z]' '[A-Z]')

#  Command line arguments.
export SIGINP=${1:-${SIGINP}}
export NEMSINP=${NEMSINP:-$SIGINP}
export FLXINP=${2:-${FLXINP}}
export FLXIOUT=${3:-${FLXIOUT}}
export PGBOUT=${4:-${PGBOUT}}
#export PGIOUT=${5:-${PGIOUT}}
export PGIOUT=${PGIOUT:-pgb.idx}
export IO=${6:-${IO:-0}}
export JO=${7:-${JO:-0}}
export IGEN=${8:-${IGEN:-0}}
export CCPOST=${CCPOST:-NO}
#  Directories.
export NWPROD=${NWPROD:-/nwprod}
export EXECUTIL=${EXECUTIL:-$NWPROD/util/exec}
export USHUTIL=${USHUTIL:-$NWPROD/util/ush}
export EXECGLOBAL=${EXECGLOBAL:-$NWPROD/exec}
export USHGLOBAL=${USHGLOBAL:-$NWPROD/ush}
export DATA=${DATA:-$(pwd)}
#  Filenames.
export MP=${MP:-$([[ $LOADL_STEP_TYPE = PARALLEL ]]&&echo "p"||echo "s")}
export XC=${XC}
export POSTGPEXEC=${POSTGPEXEC:-$EXECGLOBAL/ncep_post}
export OVERPARMEXEC=${OVERPARMEXEC:-$EXECGLOBAL/overparm_grib}
export GRBINDEX=${GRBINDEX:-$EXECUTIL/grbindex$XC}
export GRB2INDEX=${GRB2INDEX:-$EXECUTIL/grb2index$XC}
export ANOMCATSH=${ANOMCATSH:-$USHGLOBAL/global_anomcat.sh}
export CHGRESSH=${CHGRESSH:-$USHGLOBAL/global_chgres.sh}
export POSTGPLIST=${POSTGPLIST:-/dev/null}
export INISCRIPT=${INISCRIPT:-""}
export ERRSCRIPT=${ERRSCRIPT:-'eval [[ $err = 0 ]]'}
export LOGSCRIPT=${LOGSCRIPT:-""}
export ENDSCRIPT=${ENDSCRIPT:-}
export GFSOUT=${GFSOUT:-gfsout}
export PARMSUBDA=${PARMSUBDA:-parms/parm_am}
export PARMGLOBAL=${PARMGLOBAL:-$NWPROD/$PARMSUBDA}
export CTLFILE=${CTLFILE:-$PARMGLOBAL/gfs_cntrl.parm}
export MODEL_OUT_FORM=${MODEL_OUT_FORM:-grib}
export GRIBVERSION=${GRIBVERSION:-grib1}
#  Other variables.
export POSTGPVARS=$POSTGPVARS
export NTHREADS=${NTHREADS:-1}
export PGMOUT=${PGMOUT:-${pgmout:-'&1'}}
export PGMERR=${PGMERR:-${pgmerr:-'&2'}}
export CHGRESTHREAD=${CHGRESTHREAD:-1}
export FILTER=${FILTER:-1}
export GENPSICHI=${GENPSICHI:-NO}
export GENPSICHIEXE=${GENPSICHIEXE:-$EXECGLOBAL/genpsiandchi}
export BMPYXML=${BMPYXML:-FLAT}
#export D3DINP=${D3DINP:-/dev/null}

typeset -L1 l=$PGMOUT
[[ $l = '&' ]]&&a=''||a='>'
export REDOUT=${REDOUT:-'1>'$a}
typeset -L1 l=$PGMERR
[[ $l = '&' ]]&&a=''||a='>'
export REDERR=${REDERR:-'2>'$a}

#if [ $machine = IBMP6 -o $machine = WCOSS_C ] ; then
#  typeset -L1 l=$PGMOUT
#  [[ $l = '&' ]]&&a=''||a='>'
#  export REDOUT=${REDOUT:-'1>'$a}
#  typeset -L1 l=$PGMERR
#  [[ $l = '&' ]]&&a=''||a='>'
#  export REDERR=${REDERR:-'2>'$a}
#else
#  export REDOUT=${REDOUT:-'1>'}
#  export REDERR=${REDERR:-'2>'}
#fi
################################################################################
#  Preprocessing
$INISCRIPT

# Chuang: Run chgres if OUTTYP=1 or 0

export IDRT=${IDRT:-4}

if [ $OUTTYP -le 3 ] ; then
 if [ ! -s $SIGINP ] ; then      # exit if SIGINP does not exist
   echo "sigma file not found, exitting"
   exit 111
 fi
 export SIGHDR=${SIGHDR:-$NWPROD/exec/global_sighdr}
 export LONB=${LONB:-$(echo lonb|$SIGHDR $SIGINP)}
 export LATB=${LATB:-$(echo latb|$SIGHDR $SIGINP)}
else
 if [ ! -s $NEMSINP ] ; then
   echo "nemsio file not found, exitting"
   exit 112
 fi
 export SIGHDR=${SIGHDR:-$NWPROD/util/exec/nemsio_get}
 export JCAP=${LONF:-$($SIGHDR ${NEMSINP}$FM JCAP|grep -i "lonf" |awk -F"= " '{print $2}' |awk -F" " '{print $1}')}
 export LONF=${LONF:-$($SIGHDR ${NEMSINP}$FM DIMX|grep -i "dimx" |awk -F"= " '{print $2}' |awk -F" " '{print $1}')}
 export LATG=${LATG:-$($SIGHDR ${NEMSINP}$FM DIMY|grep -i "dimy" |awk -F"= " '{print $2}' |awk -F" " '{print $1}')}
 export LONB=${LONB:-LONF}
 export LATB=${LATB:-LATG}
 if [ $IDRT -eq 4 ] ; then
   export LONSPERLAT=${LONSPERLAT:-$NWPROD/fix/fix_am/global_lonsperlat.t$JCAP.$LONF.$LATG.txt}
 fi
fi

export SIGLEVEL=${SIGLEVEL:-NULL}

if [ $OUTTYP -le 1 ] ; then
 export JCAP=${JCAP:-$(echo jcap|$SIGHDR $SIGINP)}
 export LEVS=${LEVS:-$(echo levs|$SIGHDR $SIGINP)}
 export IDVC=${IDVC:-$(echo idvc|$SIGHDR $SIGINP)}
 export IDVM=${IDVM:-$(echo idvm|$SIGHDR $SIGINP)}
 export NVCOORD=${NVCOORD:-$(echo nvcoord|$SIGHDR $SIGINP)}
 export IVSSIG=${IVSSIG:-$(echo ivs|$SIGHDR $SIGINP)}
 export LATCH=${LATCH:-8}
 if [ $OUTTYP -eq 1 ] ; then 
  export CHGRESVARS="IDVC=$IDVC,IDVM=$IDVM,NVCOORD=$NVCOORD,IVSSIG=$IVSSIG,LATCH=$LATCH,"  
 elif [ $OUTTYP -eq 0 ] ; then
  export CHGRESVARS="LATCH=$LATCH,$CHGRESVARS"
 fi 
#export SIGLEVEL=${SIGLEVEL:-""}
 #export SIGLEVEL=${SIGLEVEL:-"$NWPROD/fix/global_hyblev.l${LEVS}.txt"}
 # specify threads for running chgres
 export OMP_NUM_THREADS=$CHGRESTHREAD 
 export NTHREADS=$OMP_NUM_THREADS
 export APRUNC=${APRUNCP:-${APRUNC:-""}}

 $CHGRESSH

 export ERR=$?
 export err=$ERR
 $ERRSCRIPT||exit 1

 if [ $CCPOST = YES ] ; then sleep 1 ; fi
 export MODEL_OUT_FORM=grib

elif [ $OUTTYP -eq 3 ] ; then   #run post to read sigma file directly if OUTTYP=3
 export MODEL_OUT_FORM=sigio
 export GFSOUT=$SIGINP
elif [ $OUTTYP -eq 4 ] ; then   #run post to read nemsio file if OUTTYP=4
 export MPIIO=${MPIIO:-mpiio}
 export MODEL_OUT_FORM=${MODEL_OUT_FORM:-binarynemsio$MPIIO}
 export GFSOUT=$NEMSINP
fi
#
if [ $machine = THEIA ] ; then
 export MPICH_FAST_MEMCPY=${MPICH_FAST_MEMCPY:-"ENABLE"}
 export MPI_BUFS_PER_PROC=${MPI_BUFS_PER_PROC:-2048}
 export MPI_BUFS_PER_HOST=${MPI_BUFS_PER_HOST:-2048}
#export MPI_IB_RAILS=${MPI_IB_RAILS:-2}
 export MPI_GROUP_MAX=${MPI_GROUP_MAX:-128}
 export KMP_STACKSIZE=${KMP_STACKSIZE:-1024000}
#export MPI_STATS=${MPI_STATS:-1}
#export MPI_STATS_FILE=${MPI_STATS_FILE:-$(pwd)/mpi_stats.$PBS_JOBID}
#export MPI_MEMMAP_OFF=${MPI_MEMMAP_OFF:-1}
 export OMP_NUM_THREADS=${OMP_NUM_THREADSP:-${OMP_NUM_THREADS:-8}}
 export APRUN=${APRUNP:-${APRUN_NP:-mpirun}}
elif [ $machine = WCOSS ] ; then
 export MP_EUIDEVICE=${MP_EUIDEVICE:-min}
 export KMP_STACKSIZE=${KMP_STACKSIZE:-1024000}
 export MPICH_ALLTOALL_THROTTLE=${MPICH_ALLTOALL_THROTTLE:-0}
 export MP_SINGLE_THREAD=${MP_SINGLE_THREAD:-yes}
#export MP_EUILIB=${MP_EUILIB:-us}
#export MP_EAGER_LIMIT=${MP_EAGER_LIMIT:-64K}
#export FORT_BUFFERED=${FORT_BUFFERED:-true}
 export OMP_NUM_THREADS=${OMP_NUM_THREADSP:-${OMP_NUM_THREADS:-8}}
 export APRUN=${APRUNP:-${APRUN_NP:-mpirun.lsf}}
elif [ $machine = WCOSS_C ] ; then
 export KMP_STACKSIZE=${KMP_STACKSIZE:-1024000}
 export OMP_NUM_THREADS=${OMP_NUM_THREADSP:-${OMP_NUM_THREADS:-8}}
 export OMP_NUM_THREADS_ANOM=${OMP_NUM_THREA_ANOM:-8}
 export APRUN=${APRUNP:-${APRUN_NP:-aprun}}
 export preanom=${preanom:-"aprun -n1 -N1 -d$OMP_NUM_THREADS_ANOM "}
#module load iobuf
#export IOBUF_PARAMS=${IOBUF_PARAMS:-'*:size=8M:verbose'}
#export IOBUF_PARAMS=${IOBUF_PARAMS:-'*:size=8M:%stdout:size=2M:verbose'}
elif [ $machine = GAEA ] ; then
 export KMP_STACKSIZE=${KMP_STACKSIZE:-1024000}
 export OMP_NUM_THREADS=${OMP_NUM_THREADSP:-${OMP_NUM_THREADS:-8}}
 export OMP_NUM_THREADS_ANOM=${OMP_NUM_THREA_ANOM:-8}
 export APRUN=${APRUNP:-${APRUN_NP:-aprun}}
 if [ $OMP_NUM_THREADS_ANOM -eq 1 ] ; then
  export preanom=${preanom:-"aprun -n1 -N1 "}
 else
  export preanom=${preanom:-"aprun -n1 -N1 -d$((OMP_NUM_THREADS_ANOM-1)) "}
 fi
 echo 'Not finished yet'

fi

export POST_THREADS=${POST_THREADS:-$OMP_NUM_THREADS}
export APRUN=${APRUNP:-${APRUN:-""}}
export precpgb=${precpgb:-""}
export preanom=${preanom:-$precpgb}

pwd=$(pwd)
if [[ -d $DATA ]] ; then
   mkdata=NO
else
   mkdir -p $DATA
   mkdata=YES
fi
cd $DATA||exit 99

if [ $IDRT -eq 4 ] ; then
  export LONSPERLAT=${LONSPERLAT:-${fixdir:-$NWPROD/fix/fix_am}/global_lonsperlat.t$JCAP.$LONB.$LATB.txt}
  cp $LONSPERLAT ./lonsperlat.dat
fi

################################################################################
#  Post GRIB

export PGM=$POSTGPEXEC
export pgm=$PGM
$LOGSCRIPT
cat <<EOF >postgp.inp.nml$$
 &NAMPGB
 $POSTGPVARS
/
EOF

if [[ $VERBOSE = YES ]] ; then
   cat postgp.inp.nml$$
fi

# making the time stamp format for ncep post
export YY=$(echo $VDATE | cut -c1-4)
export MM=$(echo $VDATE | cut -c5-6)
export DD=$(echo $VDATE | cut -c7-8)
export HH=$(echo $VDATE | cut -c9-10)

cat  > itag <<EOF
$GFSOUT
$MODEL_OUT_FORM
$GRIBVERSION
${YY}-${MM}-${DD}_${HH}:00:00
GFS
$FLXINP
$D3DINP
EOF

cat postgp.inp.nml$$ >> itag

cat itag

rm -f fort.*

#ln -sf $SIGINP     postgp.inp.sig$$
#ln -sf $FLXINP     postgp.inp.flx$$
#ln -sf $PGBOUT     postgp.out.pgb$$

# change model generating Grib number 

if [ $GRIBVERSION = grib1 ] ; then
  if [ $IGEN -le 9 ] ; then
   cat $CTLFILE | sed s:00082:0000${IGEN}: > ./gfs_cntrl.parm
  elif [ $IGEN -le 99 ] ; then
   cat $CTLFILE | sed s:00082:000${IGEN}:  > ./gfs_cntrl.parm
  elif [ $IGEN -le 999 ] ; then
   cat $CTLFILE | sed s:00082:00${IGEN}:   > ./gfs_cntrl.parm
  else
   ln -sf $CTLFILE ./gfs_cntrl.parm 
  fi 
  ln -sf ./gfs_cntrl.parm fort.14

elif [ $GRIBVERSION = grib2 ] ; then
  export POSTGRB2TBL=${POSTGRB2TBL:-${NWPROD:-/nwprod}/lib/sorc/g2tmpl/params_grib2_tbl_new}
  cp $POSTGRB2TBL .

  if [ $BMPYXML = 'FLAT' ] ; then
    export PostFlatFile=${PostFlatFile:-$PARMGLOBAL/postxconfig-NT-GFS.txt}
    cp $PostFlatFile ./postxconfig-NT.txt
    if [ $ens = YES ] ; then
     sed < $PostFlatFile -e "s#negatively_pert_fcst#${ens_pert_type}#" > ./postxconfig-NT.txt
    fi
  else
    export POSTAVBLFLD=${POSTAVBLFLD:-$PARMGLOBAL/post_avblflds.xml}
    cp $POSTAVBLFLD .
    cp $CTLFILE     postcntrl.xml
  fi
fi
CTL=$(basename $CTLFILE)

ln -sf griddef.out fort.110
if [ $IDRT -eq 4 ] ; then
  ln -sf $LONSPERLAT lonsperlat.dat
else
  ln -sf /dev/null lonsperlat.dat
fi
MICRO_PHYS_DATA=${MICRO_PHYS_DATA:-$NWPROD/$PARMSUBDA/nam_micro_lookup.dat}
cp $MICRO_PHYS_DATA ./eta_micro_lookup.dat

if [ $SIGLEVEL != NULL ] ; then
 cp $SIGLEVEL hyblev_file
else
 hyblev_file=/dev/null
fi

eval $APRUN $POSTGPEXEC < itag > outpost_gfs_${VDATE}_${CTL}


export ERR=$?
export err=$ERR
$ERRSCRIPT||exit 2


# Filter SLP and 500 mb height using copygb, change GRIB ID, and then
# cat the filtered fields to the pressure GRIB file, from Iredell

if [ $FILTER = "1" ] ; then

  if [ $GRIBVERSION = grib1 ] ; then
    copygb=${COPYGB:-$EXECUTIL/copygb}

#   $copygb -x -i'4,0,80' -k'4*-1,1,102' $PGBOUT tfile
    ${precpgb}$copygb -x -i'4,0,80' -k'4*-1,1,102' $PGBOUT tfile
    ln -s -f tfile fort.11
    ln -s -f prmsl fort.51
    echo 0 2|$OVERPARMEXEC

#   $copygb -x -i'4,1,5' -k'4*-1,7,100,500' $PGBOUT tfile
    ${precpgb}$copygb -x -i'4,1,5' -k'4*-1,7,100,500' $PGBOUT tfile
    ln -s -f tfile fort.11
    ln -s -f h5wav fort.51
    echo 0 222|$OVERPARMEXEC

#   cat $PGBOUT prmsl h5wav >> $PGBOUT
    cat  prmsl h5wav >> $PGBOUT

  elif [ $GRIBVERSION = grib2 ] ; then

    copygb2=${COPYGB2:-$EXECUTIL/copygb2}
    wgrib2=${WGRIB2:-$EXECUTIL/wgrib2}

    ${precpgb}$copygb2 -x -i'4,0,80' -k'0 3 0 7*-9999 101 0 0' $PGBOUT tfile
    ${precpgb}$wgrib2 tfile -set_byte 4 11 1 -grib prmsl

    ${precpgb}$copygb2 -x -i'4,1,5' -k'0 3 5 7*-9999 100 0 50000' $PGBOUT tfile
    ${precpgb}$wgrib2 tfile -set_byte 4 11 193 -grib h5wav

#  cat $PGBOUT prmsl h5wav >> $PGBOUT
    cat  prmsl h5wav >> $PGBOUT
  fi
fi

################################################################################
#  Anomaly concatenation
#  for now just do anomaly concentration for grib1

if [ $GRIBVERSION = grib1 ] ; then
 if [[ -x $ANOMCATSH ]] ; then
   if [[ -n $PGIOUT ]] ; then
      $GRBINDEX $PGBOUT $PGIOUT
   fi
   export PGM=$ANOMCATSH
   export pgm=$PGM
   $LOGSCRIPT

   eval $ANOMCATSH $PGBOUT $PGIOUT $REDOUT$PGMOUT $REDERR$PGMERR

   export ERR=$?
   export err=$ERR
   $ERRSCRIPT||exit 3
 fi
fi
################################################################################
#  Make GRIB index file
if [[ -n $PGIOUT ]] ; then
  if [ $GRIBVERSION = grib2 ] ; then
     $GRB2INDEX $PGBOUT $PGIOUT
   else
     $GRBINDEX $PGBOUT $PGIOUT
  fi
fi
if [[ -r $FLXINP && -n $FLXIOUT ]] ; then
   $GRBINDEX $FLXINP $FLXIOUT
fi
################################################################################
# generate psi and chi
echo "GENPSICHI= " $GENPSICHI
if [ $GENPSICHI = YES ] ; then
#echo "PGBOUT PGIOUT=" $PGBOUT $PGIOUT
#echo "YY MM=" $YY $MM
 export OMP_NUM_THREADS=$CHGRESTHREAD 
 export psichifile=./psichi.grb
 $GENPSICHIEXE < postgp.inp.nml$$
 rc=$?
 if [[ $rc -ne 0 ]] ; then
   echo 'Nonzero return code rc= '$rc
   if [ ${RUN_ENVIR:-dev} = dev ] ; then
      exit 3
   else
      export err=$rc; $DATA/err_chk
   fi
 fi
 cat ./psichi.grb >> $PGBOUT
fi
################################################################################
#  Postprocessing
cd $pwd
[[ $mkdata = YES ]]&&rmdir $DATA
$ENDSCRIPT
set +x
if [[ $VERBOSE = YES ]] ; then
   echo $(date) EXITING $0 with return code $err >&2
fi
exit $err
