#makefile to build io_int that does binary i/o

OBJS    = io_int.o io_int_idx.o module_io_int_idx.o pack_utils.o
LIB     = libwrfio_int.a

CPP1    = $(CPP) $(TRADFLAG) $(ARCHFLAGS)
M4      = m4 -Uinclude -Uindex -Ulen

.SUFFIXES: .F90 .F .f .f90 .c .h .o .code

all : $(LIB) 

$(LIB): $(OBJS)
	$(RM) $(LIB)
	$(AR) $(ARFLAGS) $(LIB) $(OBJS)
	$(RANLIB) $(LIB)

io_int.f: io_int.F90 module_internal_header_util.o
	$(CPP1) $*.F90 | $(M4) - > $@

io_int.o: io_int.f intio_tags.h
	$(FC) $(FCFLAGS) -o $@ -c $*.f

io_int_idx_tags.h: intio_tags.h
	awk '{print "#define", toupper($$4), $$6}' < intio_tags.h > $@

io_int_idx.o: io_int_idx.c io_int_idx.h io_int_idx_tags.h
	$(CC) -o $@ -c $(CFLAGS_LOCAL) $*.c

module_io_int_idx.o: module_io_int_idx.f
	$(FC) $(FCFLAGS) -o $@ -c $*.f

module_io_int_idx.f: module_io_int_idx.F90
	$(CPP1) $*.F90 > $@

module_internal_header_util.f: module_internal_header_util.F pack_utils.o
	-$(RM) $@
	$(CPP1) module_internal_header_util.F > $@

module_internal_header_util.o: module_internal_header_util.f
	$(FC) $(FCFLAGS) -o $@ -c $*.f

pack_utils.o :
	$(CC) -c $(CFLAGS) pack_utils.c

superclean:
	-$(RM) *.f *.o *.obj *.i *.mod $(LIB) io_int_idx_tags.h \
		io_int_idx_tags.h
