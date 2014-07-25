################################################################################
# This makefile was produced by fmjif at 17:21:11 on 97/02/24.
SHELL=	/bin/sh
CMD=	../../../exec/anomgb
FOPTS=  $(FFLAGSM)
LOPTS=  $(FFLAGS)
INCS=	
OBJS=	anomgb.o
LIBS=  $(LIBSMOD)	
################################################################################
$(CMD):	$(OBJS)
	$(FCMP) $(LOPTS) $(OBJS) $(LIBS) -o $(CMD)
$(OBJS):	$(INCS)
anomgb.o:	anomgb.f
	$(FCMP) -c $(FOPTS) $*.f
clean:	
	-rm -f $(OBJS)
