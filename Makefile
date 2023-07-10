GCLDIR := ./gclib
SEARCHDIRS := -I${GCLDIR}
# Use the correct compiler (CXX, not CC) and respect the environment
# by using ?=
CXX    ?= g++

BASEFLAGS  = -Wall ${SEARCHDIRS} $(MARCH) -DENABLE_COMPRESSION=1 -D_FILE_OFFSET_BITS=64 \
-D_LARGEFILE_SOURCE -fno-exceptions -fno-rtti -fno-strict-aliasing \
-D_REENTRANT


ifeq ($(findstring debug,$(MAKECMDGOALS)),)
  DBGFLAGS = -O2 -g -DNDEBUG
else
  DBGFLAGS = -g -DDEBUG -D_DEBUG
  LDFLAGS += -g
endif

ifeq ($(findstring nommap,$(MAKECMDGOALS)),)
  CXXFLAGS += $(DBGFLAGS) $(BASEFLAGS)
else
  CXXFLAGS += $(DBGFLAGS) $(BASEFLAGS) -DNO_MMAP
endif

%.o : %.c
	${CXX} ${CXXFLAGS} -c $< -o $@

%.o : %.cc
	${CXX} ${CXXFLAGS} -c $< -o $@

%.o : %.C
	${CXX} ${CXXFLAGS} -c $< -o $@

%.o : %.cpp
	${CXX} ${CXXFLAGS} -c $< -o $@

%.o : %.cxx
	${CXX} ${CXXFLAGS} -c $< -o $@

# C/C++ linker

LINKER    := ${CXX}
#LDFLAGS = 
#uncomment this when ENABLE_COMPRESSION
LDFLAGS    += -lz

.PHONY : all
all:    cdbfasta cdbyank 
debug:  all
nommap: all
#when compression is enabled:
#cdbfasta:  ./cdbfasta.o ./gcdbz.o  ...
cdbfasta:  ./cdbfasta.o ./gcdbz.o ${GCLDIR}/gcdb.o ${GCLDIR}/GBase.o ${GCLDIR}/GStr.o ${GCLDIR}/GArgs.o
	${LINKER} -o $@ ${filter-out %.a %.so, $^} $(LDFLAGS)
#cdbyank :  ./cdbyank.o ./gcdbz.o
cdbyank :  ./cdbyank.o ./gcdbz.o ${GCLDIR}/gcdb.o ${GCLDIR}/GBase.o ${GCLDIR}/GStr.o ${GCLDIR}/GArgs.o
	${LINKER} -o $@ ${filter-out %.a %.so, $^} $(LDFLAGS)

# target for removing all object files

.PHONY : tidy
tidy::
	@${RM} core cdbfasta cdbyank *.o ${GCLDIR}/gcdb.o ${GCLDIR}/GBase.o ${GCLDIR}/GStr.o ${GCLDIR}/GArgs.o

# target for removing all object files

.PHONY : clean
clean:: tidy
	@${RM} core cdbfasta cdbyank *.o ${GCLDIR}/gcdb.o ${GCLDIR}/GBase.o ${GCLDIR}/GStr.o ${GCLDIR}/GArgs.o


