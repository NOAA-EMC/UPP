SHELL  =  /bin/sh
LIBS   =  $(LIBSMOD)

FFLAGS =  $(FFLAGSM)

LDFLAGS=  $(FFLAGS)

SRCS=   overparm_grib.f
OBJS=   overparm_grib.o

FC=     $(FCMP)

CMD=    ../../exec/overparm_grib

all:    $(CMD)

$(CMD): $(OBJS)
	$(FC) $(LDFLAGS) -o $(@) $(OBJS) $(LIBS)
clean:
	-rm -f $(OBJS)
