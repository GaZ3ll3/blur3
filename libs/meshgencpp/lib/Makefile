### definitions ###

# default/library definitions
include $(MESHGENCPP)/lib/Makefile.defs

### targets ###

.cpp.o: ; $(LINK) $(FFLAGS) $*.cpp -o $*.o

lib: $(TOPLIB_OBJECTS_COMMON)

include $(MESHGENCPP)/lib/Makefile.targets

### DO NOT remove this line - make depends on it ###
