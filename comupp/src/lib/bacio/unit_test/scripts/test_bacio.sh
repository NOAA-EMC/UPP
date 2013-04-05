#!/bin/sh
 
crtdir=`pwd` 
libdir=${crtdir}/../../lib
libsrcdir=${crtdir}/../../src
srcdir=${crtdir}/../src
datadir=${crtdir}/../data
exedir=${crtdir}/../exe
confdir=${crtdir}/../conf
#
echo "BACIO unit test start ..."
rm $exedir/*
#
###########################################################
# detect machine:
case `hostid` in

  0xac7c006)          MACHINE_ID=ccs ;; ### cirrus1
  0xac7c012)          MACHINE_ID=ccs ;; ### cirrus2

  0xac7a006)          MACHINE_ID=ccs ;; ### stratus1
  0xac7a012)          MACHINE_ID=ccs ;; ### stratus2

  0xffffffffbcc090b2) MACHINE_ID=gaea ;; ### gaea1
  0xffffffffbcc08fb2) MACHINE_ID=gaea ;; ### gaea2
  0xffffffffbcc08eb2) MACHINE_ID=gaea ;; ### gaea3
  0xffffffffbcc08db2) MACHINE_ID=gaea ;; ### gaea4

  ae0a2500)           MACHINE_ID=zeus ;; ### zeus1
  ae0a2600)           MACHINE_ID=zeus ;; ### zeus2
  ae0a2700)           MACHINE_ID=zeus ;; ### zeus3
  ae0a2800)           MACHINE_ID=zeus ;; ### zeus4
  ae0a2900)           MACHINE_ID=zeus ;; ### zeus5
  ae0a2a00)           MACHINE_ID=zeus ;; ### zeus6
  ae0a7d13)           MACHINE_ID=zeus ;; ### zeus7
  ae0a7e13)           MACHINE_ID=zeus ;; ### zeus8

  b20a057e)           MACHINE_ID=jet ;; ### jet1
  b20a067e)           MACHINE_ID=jet ;; ### jet2
  b20a077e)           MACHINE_ID=jet ;; ### jet3
  b20a087e)           MACHINE_ID=jet ;; ### jet4

  3ac0c803)           MACHINE_ID=wcoss ;; ### current

  010afa01)           MACHINE_ID=wcoss ;; ### eddy

esac

#typeset -Z3 TEST_NR
export TEST_NR=0

echo " **************** compiling start *********************** "
if [ ${MACHINE_ID} = zeus ]; then
  echo " Compile bacio lib...."
  cd $libsrcdir
  ./makebacio_zeus.sh
  if [ -s $libdir/libbacio_4.a -a -s $libdir/libbacio_8.a ]; then
    echo " Compile bacio lib.... finished"
  else
    echo " Compile bacio lib.... failed"
    exit 8
  fi
  cd $confdir
  cp conf_zeus conf_unittest
elif [ ${MACHINE_ID} = ccs ]; then
  echo " Compile bacio lib...."
  cd $libsrcdir
  ./makebacio_ccs.sh
  if [ -s $libdir/libbacio_4.a -a -s $libdir/libbacio_8.a ]; then
    echo " Compile bacio lib.... finished"
  else
    echo " Compile bacio lib.... failed"
    exit 8
  fi
  cd $confdir
  cp conf_ccs conf_unittest
elif [ ${MACHINE_ID} = wcoss ]; then
  echo " Compile bacio lib...."
  cd $libsrcdir
  ./makebacio_wcoss.sh
  if [ -s $libdir/libbacio_4.a -a -s $libdir/libbacio_8.a ]; then
    echo " Compile bacio lib.... finished"
  else
    echo " Compile bacio lib.... failed"
    exit 8
  fi
  cd $confdir
  cp conf_current conf_unittest
fi

echo " **************** testing starts *********************** "
###########################################################
echo " "
(( TEST_NR=TEST_NR+1 ))
echo " test ${TEST_NR}: test_baciof"
cd $srcdir
make -f test_baciof_mk
if [ ! -s $exedir/test_baciof ]; then
  echo "compiling error, exit!"
  exit 8
fi
cd $exedir
./test_baciof >aa
testend=`grep "ends normally" aa |wc -l`
if [ $testend -lt 1 ]; then
  echo "test_baciof stoped!"
  exit 8
fi
cmp $datadir/data_be data_baciof_be >aa
testend=`grep "differ:" aa |wc -l`
if [ $testend -ne 0 ]; then
  echo "test_baciof stoped!"
  exit 8
fi
rm data_baciof_be aa
 

###########################################################
if [ ${MACHINE_ID} != ccs ]; then
echo " "
(( TEST_NR=TEST_NR+1 ))
echo " test ${TEST_NR}: test_baciof_le"
cd $srcdir
make -f test_baciof_le_mk
if [ ! -s $exedir/test_baciof_le ]; then
  echo "compiling error, exit!"
  exit 8
fi
cd $exedir
./test_baciof_le >aa
testend=`grep "ends normally" aa |wc -l`
if [ $testend -lt 1 ]; then
  echo "test_baciof stoped!"
  exit 8
fi
cmp $datadir/data_le data_baciof_le >aa
testend=`grep "differ:" aa |wc -l`
if [ $testend -ne 0 ]; then
  echo "test_baciof_le stoped!"
  exit 8
fi
rm data_baciof_le aa
fi

###########################################################
echo " "
(( TEST_NR=TEST_NR+1 ))
echo " test ${TEST_NR}: test_bafrof"
cd $srcdir
make -f test_bafrio_mk
if [ ! -s $exedir/test_bafrio ]; then
  echo "compiling error, exit!"
  exit 8
fi
cd $exedir
./test_bafrio >aa
testend=`grep "ends normally" aa |wc -l`
if [ $testend -lt 1 ]; then
  echo "test_bafrio stoped!"
  exit 8
fi
cmp $datadir/data_be data_bafrio_be >aa
testend=`grep "differ:" aa |wc -l`
cmp $datadir/data_be data_bafrio_be1 >bb
testend1=`grep "differ:" bb |wc -l`
testend=`expr $testend + $testend1`
if [ $testend -ne 0 ]; then
  echo "test_bafriof stoped!"
  exit 8
fi
rm data_bafrio_be data_bafrio_be1 aa bb
 

###########################################################
if [ ${MACHINE_ID} != ccs ]; then
echo " "
(( TEST_NR=TEST_NR+1 ))
echo " test ${TEST_NR}: test_bafriof_le"
cd $srcdir
make -f test_bafrio_le_mk
if [ ! -s $exedir/test_bafrio_le ]; then
  echo "compiling error, exit!"
  exit 8
fi
cd $exedir
./test_bafrio_le >aa
testend=`grep "ends normally" aa |wc -l`
if [ $testend -lt 1 ]; then
  echo "test_bafriof stoped!"
  exit 8
fi
cmp $datadir/data_le data_bafrio_le >aa
testend=`grep "differ:" aa |wc -l`
cmp $datadir/data_le data_bafrio_le1 >bb
testend1=`grep "differ:" bb |wc -l`
testend=`expr $testend + $testend1`
if [ $testend -ne 0 ]; then
  echo "test_bafriof_le stoped!"
  exit 8
fi
rm data_bafrio_le data_bafrio_le1 aa bb
fi

###########################################################
echo " "
(( TEST_NR=TEST_NR+1 ))
echo " test ${TEST_NR}: test_wryte"
cd $srcdir
make -f test_wryte_mk
if [ ! -s $exedir/test_wryte ]; then
  echo "compiling error, exit!"
  exit 8
fi
cd $exedir
./test_wryte >aa
testend=`grep "ends normally" aa |wc -l`
if [ $testend -lt 1 ]; then
  echo "test_bafriof stoped!"
  exit 8
fi
cmp $datadir/data_wryte data_wryte >aa
testend=`grep "differ:" aa |wc -l`
if [ $testend -ne 0 ]; then
  echo "test_bafriof stoped!"
  exit 8
fi
rm data_wryte aa

rm $exedir/*
###########################################################
echo " "
echo "BACIO lib unit test finished!"
