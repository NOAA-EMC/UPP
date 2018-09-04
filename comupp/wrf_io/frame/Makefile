#

LN      =       ln -sf
MAKE    =       make -i -r
RM      =       rm -f



OBJS    =       \
                wrf_num_bytes_between.o    \
                libmassv.o                 \
                collect_on_comm.o          \
                hires_timer.o

#compile as a .o but do not link into the main library
SPECIAL =       pack_utils.o
   

include ../wrf_makerules

TARGETDIR    =  ./

wrf_num_bytes_between.o :
	$(CC) -c $(CFLAGS) wrf_num_bytes_between.c

pack_utils.o :
	$(CC) -c $(CFLAGS) -DIWORDSIZE=$(IWORDSIZE) pack_utils.c

clean:
	@ echo 'use the clean script'

md_calls.inc : md_calls.m4
	if [ "$(M4)" = NA ] ; then \
	  /bin/cp ../arch/md_calls.inc . ; \
	else \
	  $(M4) md_calls.m4 > md_calls.inc ; \
	fi

# DEPENDENCIES : only dependencies after this line (don't remove the word DEPENDENCIES)

# DO NOT DELETE
