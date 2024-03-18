#!/bin/bash

case $(hostname -f) in

  alogin01.acorn.wcoss2.ncep.noaa.gov)  MACHINE_ID=wcoss2_a ;; ### acorn
  alogin02.acorn.wcoss2.ncep.noaa.gov)  MACHINE_ID=wcoss2_a ;; ### acorn
  adecflow01.acorn.wcoss2.ncep.noaa.gov) MACHINE_ID=wcoss2 ;; ### acorn
  adecflow02.acorn.wcoss2.ncep.noaa.gov) MACHINE_ID=wcoss2 ;; ### acorn
  dlogin01.dogwood.wcoss2.ncep.noaa.gov) MACHINE_ID=wcoss2 ;; ### dodwood
  dlogin02.dogwood.wcoss2.ncep.noaa.gov) MACHINE_ID=wcoss2 ;; ### dodwood
  dlogin03.dogwood.wcoss2.ncep.noaa.gov) MACHINE_ID=wcoss2 ;; ### dodwood
  dlogin04.dogwood.wcoss2.ncep.noaa.gov) MACHINE_ID=wcoss2 ;; ### dodwood
  dlogin05.dogwood.wcoss2.ncep.noaa.gov) MACHINE_ID=wcoss2 ;; ### dodwood
  dlogin06.dogwood.wcoss2.ncep.noaa.gov) MACHINE_ID=wcoss2 ;; ### dodwood
  dlogin07.dogwood.wcoss2.ncep.noaa.gov) MACHINE_ID=wcoss2 ;; ### dodwood
  dlogin08.dogwood.wcoss2.ncep.noaa.gov) MACHINE_ID=wcoss2 ;; ### dodwood
  dlogin09.dogwood.wcoss2.ncep.noaa.gov) MACHINE_ID=wcoss2 ;; ### dodwood
  clogin01.cactus.wcoss2.ncep.noaa.gov) MACHINE_ID=wcoss2 ;; ### cactus
  clogin02.cactus.wcoss2.ncep.noaa.gov) MACHINE_ID=wcoss2 ;; ### cactus
  clogin03.cactus.wcoss2.ncep.noaa.gov) MACHINE_ID=wcoss2 ;; ### cactus
  clogin04.cactus.wcoss2.ncep.noaa.gov) MACHINE_ID=wcoss2 ;; ### cactus
  clogin05.cactus.wcoss2.ncep.noaa.gov) MACHINE_ID=wcoss2 ;; ### cactus
  clogin06.cactus.wcoss2.ncep.noaa.gov) MACHINE_ID=wcoss2 ;; ### cactus
  clogin07.cactus.wcoss2.ncep.noaa.gov) MACHINE_ID=wcoss2 ;; ### cactus
  clogin08.cactus.wcoss2.ncep.noaa.gov) MACHINE_ID=wcoss2 ;; ### cactus
  clogin09.cactus.wcoss2.ncep.noaa.gov) MACHINE_ID=wcoss2 ;; ### cactus

  gaea5[1-8]*)                   MACHINE_ID=gaea ;; ### gaea9
  gaea5[1-8].ncrc.gov)           MACHINE_ID=gaea ;; ### gaea9

  hfe01)                   MACHINE_ID=hera_c ;; ### hera01
  hfe02)                   MACHINE_ID=hera_c ;; ### hera02
  hfe03)                   MACHINE_ID=hera_c ;; ### hera03
  hfe04)                   MACHINE_ID=hera_c ;; ### hera04
  hfe05)                   MACHINE_ID=hera_c ;; ### hera05
  hfe06)                   MACHINE_ID=hera_c ;; ### hera06
  hfe07)                   MACHINE_ID=hera_c ;; ### hera07
  hfe08)                   MACHINE_ID=hera_c ;; ### hera08
  hfe09)                   MACHINE_ID=hera ;; ### hera09
  hfe10)                   MACHINE_ID=hera ;; ### hera10
  hfe11)                   MACHINE_ID=hera ;; ### hera11
  hfe12)                   MACHINE_ID=hera ;; ### hera12

  fe1)                     MACHINE_ID=jet_c ;; ### jet01
  fe2)                     MACHINE_ID=jet_c ;; ### jet02
  fe3)                     MACHINE_ID=jet_c ;; ### jet03
  fe4)                     MACHINE_ID=jet_c ;; ### jet04
  fe5)                     MACHINE_ID=jet ;; ### jet05
  fe6)                     MACHINE_ID=jet ;; ### jet06
  fe7)                     MACHINE_ID=jet ;; ### jet07
  fe8)                     MACHINE_ID=jet ;; ### jet08
  tfe1)                    MACHINE_ID=jet ;; ### jet09
  tfe2)                    MACHINE_ID=jet ;; ### jet10

  Orion-login-1.HPC.MsState.Edu) MACHINE_ID=orion ;; ### orion1
  Orion-login-2.HPC.MsState.Edu) MACHINE_ID=orion ;; ### orion2
  Orion-login-3.HPC.MsState.Edu) MACHINE_ID=orion ;; ### orion3
  Orion-login-4.HPC.MsState.Edu) MACHINE_ID=orion ;; ### orion4

  Hercules-login-1.HPC.MsState.Edu) MACHINE_ID=hercules ;; ### hercules1
  Hercules-login-2.HPC.MsState.Edu) MACHINE_ID=hercules ;; ### hercules2
  Hercules-login-3.HPC.MsState.Edu) MACHINE_ID=hercules ;; ### hercules3
  Hercules-login-4.HPC.MsState.Edu) MACHINE_ID=hercules ;; ### hercules4

  cheyenne1.cheyenne.ucar.edu) MACHINE_ID=cheyenne ;; ### cheyenne1
  cheyenne2.cheyenne.ucar.edu) MACHINE_ID=cheyenne ;; ### cheyenne2
  cheyenne3.cheyenne.ucar.edu) MACHINE_ID=cheyenne ;; ### cheyenne3
  cheyenne4.cheyenne.ucar.edu) MACHINE_ID=cheyenne ;; ### cheyenne4
  cheyenne5.cheyenne.ucar.edu) MACHINE_ID=cheyenne ;; ### cheyenne5
  cheyenne6.cheyenne.ucar.edu) MACHINE_ID=cheyenne ;; ### cheyenne6
  cheyenne1.ib0.cheyenne.ucar.edu) MACHINE_ID=cheyenne ;; ### cheyenne1
  cheyenne2.ib0.cheyenne.ucar.edu) MACHINE_ID=cheyenne ;; ### cheyenne2
  cheyenne3.ib0.cheyenne.ucar.edu) MACHINE_ID=cheyenne ;; ### cheyenne3
  cheyenne4.ib0.cheyenne.ucar.edu) MACHINE_ID=cheyenne ;; ### cheyenne4
  cheyenne5.ib0.cheyenne.ucar.edu) MACHINE_ID=cheyenne ;; ### cheyenne5
  cheyenne6.ib0.cheyenne.ucar.edu) MACHINE_ID=cheyenne ;; ### cheyenne6

  login1.stampede2.tacc.utexas.edu) MACHINE_ID=stampede ;; ### stampede1
  login2.stampede2.tacc.utexas.edu) MACHINE_ID=stampede ;; ### stampede2
  login3.stampede2.tacc.utexas.edu) MACHINE_ID=stampede ;; ### stampede3
  login4.stampede2.tacc.utexas.edu) MACHINE_ID=stampede ;; ### stampede4

  s4-submit.ssec.wisc.edu) MACHINE_ID=s4 ;; ### S4

  *) MACHINE_ID=unknown
esac

# Overwrite auto-detect with RT_MACHINE if set
MACHINE_ID=${RT_MACHINE:-${MACHINE_ID}}

# Append compiler
#if [ $MACHINE_ID = orion ] || [ $MACHINE_ID = hera ] || [ $MACHINE_ID = cheyenne ] || [ $MACHINE_ID = jet ] || [ $MACHINE_ID = gaea ] || [ $MACHINE_ID = stampede ] || [ $MACHINE_ID = s4 ]; then
#    MACHINE_ID=${MACHINE_ID}.${RT_COMPILER}
#fi

echo "Machine: " $MACHINE_ID 
