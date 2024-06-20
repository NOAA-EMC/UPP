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
  postxconfig-NT-GFS.txt \
  postxconfig-NT-GFS-GOES.txt \
  postxconfig-NT-GFS-ANL.txt \
  postxconfig-NT-GFS-F00.txt \
  postxconfig-NT-GFS-FLUX.txt \
  postxconfig-NT-GFS-FLUX-F00.txt \
  postxconfig-NT-NMM.txt \
  postxconfig-NT-NGAC.txt \
  postxconfig-NT-GEFS.txt \
  postxconfig-NT-GEFS-F00.txt \
  postxconfig-NT-hafs_ar_nosat.txt \
  postxconfig-NT-hafs_ar_nosat_l65.txt \
  postxconfig-NT-hafs_sat.txt \
  postxconfig-NT-hafs_nosat.txt \
  postxconfig-NT-hafs.txt \
  postxconfig-NT-GFS-TWO.txt \
  postxconfig-NT-GFS-F00-TWO.txt \
  postxconfig-NT-hrrr.txt \
  postxconfig-NT-rap.txt \
  postxconfig-NT-rrfs.txt \
  postxconfig-NT-rrfs_subh.txt \
  postxconfig-NT-UFS-aerosol.txt \
  postxconfig-NT-UFS-aerosol-F00.txt \
  postxconfig-NT-GEFS-aerosol.txt \
  postxconfig-NT-GEFS-F00-aerosol.txt \
  postxconfig-NT-AQM.txt \
  postxconfig-NT-SFS.txt

.PHONY: all clean
all: $(ALL_TXT)
clean:
	$(RM) $(ALL_TXT)

postxconfig-NT-GFS.txt:              postcntrl_gfs.xml              post_avblflds.xml
postxconfig-NT-GFS-GOES.txt:         postcntrl_gfs_goes.xml         post_avblflds.xml
postxconfig-NT-GFS-ANL.txt:          postcntrl_gfs_anl.xml          post_avblflds.xml
postxconfig-NT-GFS-F00.txt:          postcntrl_gfs_f00.xml          post_avblflds.xml
postxconfig-NT-GFS-FLUX.txt:         postcntrl_gfs_flux.xml         post_avblflds.xml
postxconfig-NT-GFS-FLUX-F00.txt:     postcntrl_gfs_flux_f00.xml     post_avblflds.xml
postxconfig-NT-NMM.txt:              nam_cntrl_cmaq.xml             nam_post_avblflds.xml
postxconfig-NT-NGAC.txt:             ngac_postcntrl.xml             ngac_post_avblflds.xml
postxconfig-NT-GEFS.txt:             postcntrl_gefs.xml             post_avblflds.xml
postxconfig-NT-GEFS-F00.txt:         postcntrl_gefs_f00.xml         post_avblflds.xml
postxconfig-NT-hafs_sat.txt:         postcntrl_hafs_sat.xml         post_avblflds.xml
postxconfig-NT-hafs_nosat.txt:       postcntrl_hafs_nosat.xml       post_avblflds.xml
postxconfig-NT-hafs_ar_nosat.txt:    postcntrl_hafs_ar_nosat.xml    post_avblflds.xml
postxconfig-NT-hafs_ar_nosat_l65.txt: postcntrl_hafs_ar_nosat_l65.xml post_avblflds.xml
postxconfig-NT-hafs.txt:             postcntrl_hafs.xml             post_avblflds.xml
postxconfig-NT-GFS-TWO.txt:          postcntrl_gfs_two.xml          post_avblflds.xml
postxconfig-NT-GFS-F00-TWO.txt:      postcntrl_gfs_f00_two.xml      post_avblflds.xml
postxconfig-NT-hrrr.txt:             hrrr_postcntrl.xml             post_avblflds_raphrrr.xml
postxconfig-NT-rap.txt:              rap_postcntrl.xml              post_avblflds_raphrrr.xml
postxconfig-NT-rrfs.txt:             rrfs_postcntrl.xml             post_avblflds.xml
postxconfig-NT-rrfs_subh.txt:        rrfs_postcntrl_subh.xml        post_avblflds.xml
postxconfig-NT-UFS-aerosol.txt:      postcntrl_ufs_aerosol.xml      post_avblflds.xml
postxconfig-NT-UFS-aerosol-F00.txt:  postcntrl_ufs_aerosol_f00.xml  post_avblflds.xml
postxconfig-NT-GEFS-aerosol.txt:     postcntrl_gefs_aerosol.xml     post_avblflds.xml
postxconfig-NT-GEFS-F00-aerosol.txt: postcntrl_gefs_aerosol_f00.xml post_avblflds.xml
postxconfig-NT-AQM.txt:              aqm.xml                        post_avblflds.xml
postxconfig-NT-SFS.txt:		     postcntrl_sfs.xml		    post_avblflds.xml

%.txt:
	$(PERL) $(PERLXML) $^ $@
