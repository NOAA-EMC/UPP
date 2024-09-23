# This makefile is to create Post flat file from available and control XML file
# Currently support GFS, GFS GOES, GFS ANL, NMM, and NGAC
# See documentation EMC POST Performance Improvement Upgrade for information
# 2015/04       Lin Gan : First version
# 2015/10       Lin Gan : Add HRRR and RAP
# 2016/03       Lin Gan : Add GEFS and new control section tag using GEFS ctrl and GFS avail file
# 2023/09       Sam Trahan : Rewritten with simpler syntax and better readability

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
  ./gfs/postxconfig-NT-gfs.txt \
  ./gfs/postxconfig-NT-gfs-goes.txt \
  ./gfs/postxconfig-NT-gfs-anl.txt \
  ./gfs/postxconfig-NT-gfs-f00.txt \
  ./gfs/postxconfig-NT-gfs-flux.txt \
  ./gfs/postxconfig-NT-gfs-flux-f00.txt \
  ./gfs/postxconfig-NT-gfs-two.txt \
  ./gfs/postxconfig-NT-gfs-f00-two.txt \
  ./gfs/postxconfig-NT-gfs-wafs-anl.txt \
  ./gfs/postxconfig-NT-gfs-wafs-ext.txt \
  ./gfs/postxconfig-NT-gfs-wafs.txt \
  ./gefs/postxconfig-NT-gefs.txt \
  ./gefs/postxconfig-NT-gefs-f00.txt \
  ./gefs/postxconfig-NT-gefs-aerosol.txt \
  ./gefs/postxconfig-NT-gefs-f00-aerosol.txt \
  ./gefs/postxconfig-NT-gefs-wafs.txt \
  ./sfs/postxconfig-NT-sfs.txt \
  postxconfig-NT-NMM.txt \
  postxconfig-NT-NGAC.txt \
  postxconfig-NT-hafs_ar_nosat.txt \
  postxconfig-NT-hafs_ar_nosat_l65.txt \
  postxconfig-NT-hafs_sat.txt \
  postxconfig-NT-hafs_nosat.txt \
  postxconfig-NT-hafs.txt \
  postxconfig-NT-hrrr.txt \
  postxconfig-NT-rap.txt \
  postxconfig-NT-rrfs.txt \
  postxconfig-NT-rrfs_subh.txt \
  postxconfig-NT-UFS-aerosol.txt \
  postxconfig-NT-UFS-aerosol-F00.txt \
  postxconfig-NT-AQM.txt

.PHONY: all clean
all: $(ALL_TXT)
clean:
	$(RM) $(ALL_TXT)

./gfs/postxconfig-NT-gfs.txt:               ./gfs/postcntrl_gfs.xml               post_avblflds.xml
./gfs/postxconfig-NT-gfs-goes.txt:          ./gfs/postcntrl_gfs_goes.xml          post_avblflds.xml
./gfs/postxconfig-NT-gfs-anl.txt:           ./gfs/postcntrl_gfs_anl.xml           post_avblflds.xml
./gfs/postxconfig-NT-gfs-f00.txt:           ./gfs/postcntrl_gfs_f00.xml           post_avblflds.xml
./gfs/postxconfig-NT-gfs-flux.txt:          ./gfs/postcntrl_gfs_flux.xml          post_avblflds.xml
./gfs/postxconfig-NT-gfs-flux-f00.txt:      ./gfs/postcntrl_gfs_flux_f00.xml      post_avblflds.xml
./gfs/postxconfig-NT-gfs-two.txt:           ./gfs/postcntrl_gfs_two.xml           post_avblflds.xml
./gfs/postxconfig-NT-gfs-f00-two.txt:       ./gfs/postcntrl_gfs_f00_two.xml       post_avblflds.xml
./gfs/postxconfig-NT-gfs-wafs-anl.txt:      ./gfs/postcntrl_gfs_wafs_anl.xml      post_avblflds.xml
./gfs/postxconfig-NT-gfs-wafs-ext.txt:      ./gfs/postcntrl_gfs_wafs_ext.xml      post_avblflds.xml
./gfs/postxconfig-NT-gfs-wafs.txt:          ./gfs/postcntrl_gfs_wafs.xml          post_avblflds.xml
./gefs/postxconfig-NT-gefs.txt:             ./gefs/postcntrl_gefs.xml             post_avblflds.xml
./gefs/postxconfig-NT-gefs-f00.txt:         ./gefs/postcntrl_gefs_f00.xml         post_avblflds.xml
./gefs/postxconfig-NT-gefs-aerosol.txt:     ./gefs/postcntrl_gefs_aerosol.xml     post_avblflds.xml
./gefs/postxconfig-NT-gefs-f00-aerosol.txt: ./gefs/postcntrl_gefs_aerosol_f00.xml post_avblflds.xml
./gefs/postxconfig-NT-gefs-wafs.txt:        ./gefs/postcntrl_gefs_wafs.xml        post_avblflds.xml
./sfs/postxconfig-NT-sfs.txt:		    ./sfs/postcntrl_sfs.xml	          post_avblflds.xml
postxconfig-NT-NMM.txt:              nam_cntrl_cmaq.xml             nam_post_avblflds.xml
postxconfig-NT-NGAC.txt:             ngac_postcntrl.xml             ngac_post_avblflds.xml
postxconfig-NT-hafs_sat.txt:         postcntrl_hafs_sat.xml         post_avblflds.xml
postxconfig-NT-hafs_nosat.txt:       postcntrl_hafs_nosat.xml       post_avblflds.xml
postxconfig-NT-hafs_ar_nosat.txt:    postcntrl_hafs_ar_nosat.xml    post_avblflds.xml
postxconfig-NT-hafs_ar_nosat_l65.txt: postcntrl_hafs_ar_nosat_l65.xml post_avblflds.xml
postxconfig-NT-hafs.txt:             postcntrl_hafs.xml             post_avblflds.xml
postxconfig-NT-hrrr.txt:             hrrr_postcntrl.xml             post_avblflds_raphrrr.xml
postxconfig-NT-rap.txt:              rap_postcntrl.xml              post_avblflds_raphrrr.xml
postxconfig-NT-rrfs.txt:             rrfs_postcntrl.xml             post_avblflds.xml
postxconfig-NT-rrfs_subh.txt:        rrfs_postcntrl_subh.xml        post_avblflds.xml
postxconfig-NT-UFS-aerosol.txt:      postcntrl_ufs_aerosol.xml      post_avblflds.xml
postxconfig-NT-UFS-aerosol-F00.txt:  postcntrl_ufs_aerosol_f00.xml  post_avblflds.xml
postxconfig-NT-AQM.txt:              aqm.xml                        post_avblflds.xml

%.txt:
	$(PERL) $(PERLXML) $^ $@
