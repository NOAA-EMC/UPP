#

LN      =       ln -sf
MAKE    =       make -i -r
RM      =       rm -f



MODULES =       module_driver_constants.o  \
		$(LLIST) \
                module_machine.o           \
                module_wrf_error.o         \
                module_sm.o                \
                module_comm_dm.o           \
                module_comm_dm_0.o         \
                module_comm_dm_1.o         \
                module_comm_dm_2.o         \
                module_comm_dm_3.o         \
                module_comm_dm_4.o         \
                module_comm_nesting_dm.o

OBJS    =       \
                wrf_num_bytes_between.o    \
                wrf_debug.o                \
                libmassv.o                 \
                collect_on_comm.o          \
                hires_timer.o

#compile as a .o but do not link into the main library
SPECIAL =       module_internal_header_util.o pack_utils.o
   

include ../wrf_makerules

TARGETDIR    =  ./

wrf_num_bytes_between.o :
	$(CC) -c $(CFLAGS) wrf_num_bytes_between.c

pack_utils.o :
	$(CC) -c $(CFLAGS) -DIWORDSIZE=$(IWORDSIZE) pack_utils.c

module_internal_header_util.o :
	$(CPP) $(CPPFLAGS) -I. module_internal_header_util.F > module_internal_header_util.f90
	$(SFC) -c $(PROMOTION) $(FCSUFFIX) $(FCDEBUG) $(FCBASEOPTS) module_internal_header_util.f90

clean:
	@ echo 'use the clean script'

md_calls.inc : md_calls.m4
	if [ "$(M4)" = NA ] ; then \
	  /bin/cp ../arch/md_calls.inc . ; \
	else \
	  $(M4) md_calls.m4 > md_calls.inc ; \
	fi

# DEPENDENCIES : only dependencies after this line (don't remove the word DEPENDENCIES)

include depend.common

# DO NOT DELETE
