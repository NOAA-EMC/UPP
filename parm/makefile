# This makefile is to create Post flat file from available and control XML file
# Currently support GFS, GFS GOES, GFS ANL, NMM, and NGAC
# See documentation EMC POST Performance Improvement Upgrade for information
# 2015/04       Lin Gan : First version
# 2015/10       Lin Gan : Add HRRR and RAP
# 2016/03       Lin Gan : Add GEFS and new control section tag using GEFS ctrl and GFS avail file
# 2023/09       Sam Trahan : Rewritten with simpler syntax and better readability
# 2025/05       Jaymes Kenyon : Adding files for MPAS applications
# 2025/05       Ben Blake : Updating files for RRFSv1
# 2025/07       Wen Meng : Add GCAFS

SHELL=/bin/sh
PERL=/usr/bin/perl
PERLXML=PostXMLPreprocessor.pl
RM=rm -f

# Rule Syntax:
#  target.txt: control.xml avblflds.xml
# Will call:
#  perl PostXMLPreprocessor.pl control.xml avblflds.xml target.txt
# Make sure all of the .txt files are in the ALL_TXT variable.

ALL_TXT= \
  ./gfs/postxconfig-NT-gfs-goes.txt \
  ./gfs/postxconfig-NT-gfs-anl.txt \
  ./gfs/postxconfig-NT-gfs-two.txt \
  ./gfs/postxconfig-NT-gfs-f00-two.txt \
  ./gefs/postxconfig-NT-gefs.txt \
  ./gefs/postxconfig-NT-gefs-f00.txt \
  ./gefs/postxconfig-NT-gefs-aerosol.txt \
  ./gefs/postxconfig-NT-gefs-f00-aerosol.txt \
  ./wafs/postxconfig-NT-gfs-wafs-anl.txt \
  ./wafs/postxconfig-NT-gfs-wafs-ext.txt \
  ./wafs/postxconfig-NT-gfs-wafs.txt \
  ./wafs/postxconfig-NT-gefs-wafs.txt \
  ./sfs/postxconfig-NT-sfs.txt \
  postxconfig-NT-NMM.txt \
  postxconfig-NT-NGAC.txt \
  postxconfig-NT-hafs_ar_nosat.txt \
  postxconfig-NT-hafs_ar_nosat_l65.txt \
  postxconfig-NT-hafs_sat.txt \
  postxconfig-NT-hafs_nosat.txt \
  postxconfig-NT-hafs.txt \
  postxconfig-NT-hrrr.txt \
  postxconfig-NT-ifi.txt \
  postxconfig-NT-rap.txt \
  ./rrfs/postxconfig-NT-rrfs.txt \
  ./rrfs/postxconfig-NT-rrfs_subh.txt \
  ./rrfs/postxconfig-NT-rrfs_f00.txt \
  ./rrfs/postxconfig-NT-rrfs_f01.txt \
  ./rrfs/postxconfig-NT-refs.txt \
  ./rrfs/postxconfig-NT-refs_f00.txt \
  ./rrfs/postxconfig-NT-refs_f01.txt \
  ./rrfs/postxconfig-NT-rrfs_firewx.txt \
  ./rrfs/postxconfig-NT-rrfs_firewx_f00.txt \
  ./rrfs/postxconfig-NT-rrfs_firewx_f01.txt \
  ./mpas/postxconfig-NT-rrfs_mpas.txt \
  ./mpas/postxconfig-NT-hfip_mpas.txt \
  ./mpas/postxconfig-NT-global_mpas.txt \
  postxconfig-NT-UFS-aerosol.txt \
  postxconfig-NT-UFS-aerosol-F00.txt \
  postxconfig-NT-AQM.txt \
  ./dafs/postxconfig-NT-hrrr_dafs.txt \
  postxconfig-NT-wps.txt \
  ./gcafs/postxconfig-NT-gcafs.txt \
  ./3drtma/postxconfig-NT-3drtma.txt

.PHONY: all clean
all: $(ALL_TXT)
clean:
	$(RM) $(ALL_TXT)

./gfs/postxconfig-NT-gfs-goes.txt:          ./gfs/postcntrl_gfs_goes.xml          post_avblflds.xml
./gfs/postxconfig-NT-gfs-anl.txt:           ./gfs/postcntrl_gfs_anl.xml           post_avblflds.xml
./gfs/postxconfig-NT-gfs-two.txt:           ./gfs/postcntrl_gfs_two.xml           post_avblflds.xml
./gfs/postxconfig-NT-gfs-f00-two.txt:       ./gfs/postcntrl_gfs_f00_two.xml       post_avblflds.xml
./gefs/postxconfig-NT-gefs.txt:             ./gefs/postcntrl_gefs.xml             post_avblflds.xml
./gefs/postxconfig-NT-gefs-f00.txt:         ./gefs/postcntrl_gefs_f00.xml         post_avblflds.xml
./gefs/postxconfig-NT-gefs-aerosol.txt:     ./gefs/postcntrl_gefs_aerosol.xml     post_avblflds.xml
./gefs/postxconfig-NT-gefs-f00-aerosol.txt: ./gefs/postcntrl_gefs_aerosol_f00.xml post_avblflds.xml
./wafs/postxconfig-NT-gfs-wafs-anl.txt:      ./wafs/postcntrl_gfs_wafs_anl.xml      post_avblflds.xml
./wafs/postxconfig-NT-gfs-wafs-ext.txt:      ./wafs/postcntrl_gfs_wafs_ext.xml      post_avblflds.xml
./wafs/postxconfig-NT-gfs-wafs.txt:          ./wafs/postcntrl_gfs_wafs.xml          post_avblflds.xml
./wafs/postxconfig-NT-gefs-wafs.txt:        ./wafs/postcntrl_gefs_wafs.xml        post_avblflds.xml
./sfs/postxconfig-NT-sfs.txt:		    ./sfs/postcntrl_sfs.xml	          post_avblflds.xml
postxconfig-NT-NMM.txt:              nam_cntrl_cmaq.xml             nam_post_avblflds.xml
postxconfig-NT-NGAC.txt:             ngac_postcntrl.xml             ngac_post_avblflds.xml
postxconfig-NT-hafs_sat.txt:         postcntrl_hafs_sat.xml         post_avblflds.xml
postxconfig-NT-hafs_nosat.txt:       postcntrl_hafs_nosat.xml       post_avblflds.xml
postxconfig-NT-hafs_ar_nosat.txt:    postcntrl_hafs_ar_nosat.xml    post_avblflds.xml
postxconfig-NT-hafs_ar_nosat_l65.txt: postcntrl_hafs_ar_nosat_l65.xml post_avblflds.xml
postxconfig-NT-hafs.txt:             postcntrl_hafs.xml             post_avblflds.xml
postxconfig-NT-hrrr.txt:             hrrr_postcntrl.xml             post_avblflds_raphrrr.xml
postxconfig-NT-ifi.txt:              ifi_postcntrl.xml              post_avblflds.xml
postxconfig-NT-rap.txt:              rap_postcntrl.xml              post_avblflds_raphrrr.xml
./rrfs/postxconfig-NT-rrfs.txt:      ./rrfs/rrfs_postcntrl.xml      post_avblflds.xml
./rrfs/postxconfig-NT-rrfs_subh.txt: ./rrfs/rrfs_postcntrl_subh.xml post_avblflds.xml
./rrfs/postxconfig-NT-rrfs_f00.txt:  ./rrfs/rrfs_postcntrl_f00.xml  post_avblflds.xml
./rrfs/postxconfig-NT-rrfs_f01.txt:  ./rrfs/rrfs_postcntrl_f01.xml  post_avblflds.xml
./rrfs/postxconfig-NT-refs.txt:      ./rrfs/refs_postcntrl.xml      post_avblflds.xml
./rrfs/postxconfig-NT-refs_f00.txt:  ./rrfs/refs_postcntrl_f00.xml  post_avblflds.xml
./rrfs/postxconfig-NT-refs_f01.txt:  ./rrfs/refs_postcntrl_f01.xml  post_avblflds.xml
./rrfs/postxconfig-NT-rrfs_firewx.txt:     ./rrfs/rrfs_postcntrl_firewx.xml     post_avblflds.xml
./rrfs/postxconfig-NT-rrfs_firewx_f00.txt: ./rrfs/rrfs_postcntrl_firewx_f00.xml post_avblflds.xml
./rrfs/postxconfig-NT-rrfs_firewx_f01.txt: ./rrfs/rrfs_postcntrl_firewx_f01.xml post_avblflds.xml
./mpas/postxconfig-NT-rrfs_mpas.txt:   ./mpas/rrfs_mpas_postcntrl.xml   post_avblflds.xml
./mpas/postxconfig-NT-hfip_mpas.txt:   ./mpas/hfip_mpas_postcntrl.xml   post_avblflds.xml
./mpas/postxconfig-NT-global_mpas.txt: ./mpas/global_mpas_postcntrl.xml post_avblflds.xml
postxconfig-NT-UFS-aerosol.txt:      postcntrl_ufs_aerosol.xml      post_avblflds.xml
postxconfig-NT-UFS-aerosol-F00.txt:  postcntrl_ufs_aerosol_f00.xml  post_avblflds.xml
postxconfig-NT-AQM.txt:              aqm.xml                        post_avblflds.xml
./dafs/postxconfig-NT-hrrr_dafs.txt: ./dafs/hrrr_postcntrl_dafs.xml post_avblflds.xml
postxconfig-NT-wps.txt:        postcntrl_wps.xml        post_avblflds.xml
./gcafs/postxconfig-NT-gcafs.txt: ./gcafs/postcntrl_gcafs.xml post_avblflds.xml
./3drtma/postxconfig-NT-3drtma.txt: ./3drtma/3drtma_postcntrl.xml post_avblflds.xml

%.txt:
	$(PERL) $(PERLXML) $^ $@
