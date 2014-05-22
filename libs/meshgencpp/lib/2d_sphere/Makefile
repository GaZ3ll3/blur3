#
#               Makefile for the lib code:
#
#       To make an executable, type:            make
#       (from the application sub-directory)

### definitions ###

# default/library definitions
include $(MESHGENCPP)/lib/2d_sphere/Makefile.defs

# compiler options

ALL_SOURCES = $(LIBSOURCES_COMMON)
ALL_OBJECTS = $(LIBOBJECTS_COMMON)

### targets ###

.cpp.o: ; $(LINK) $(FFLAGS) $*.cpp -o $*.o

lib: $(LIBOBJECTS_COMMON)

include $(MESHGENCPP)/lib/2d_sphere/Makefile.targets

cleanlib:
	rm -f *.o $(MESHGENCPP)/lib/*.o

### DO NOT remove this line - make depends on it ###