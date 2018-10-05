GCLDIR := ./gclib
SEARCHDIRS := -I${GCLDIR}
CC      := g++

BASEFLAGS  = -Wall ${SEARCHDIRS} $(MARCH) -DENABLE_COMPRESSION=1 -D_FILE_OFFSET_BITS=64 \
-D_LARGEFILE_SOURCE -fno-exceptions -fno-rtti -fno-strict-aliasing \
-D_REENTRANT


ifeq ($(findstring debug,$(MAKECMDGOALS)),)
  DBGFLAGS = -O2 -g -DNDEBUG
  LDFLAGS =
else
  DBGFLAGS = -g -DDEBUG -D_DEBUG
  LDFLAGS = -g
endif

ifeq ($(findstring nommap,$(MAKECMDGOALS)),)
  CFLAGS = $(DBGFLAGS) $(BASEFLAGS)
else
  CFLAGS = $(DBGFLAGS) $(BASEFLAGS) -DNO_MMAP
endif

%.o : %.c
	${CC} ${CFLAGS} -c $< -o $@

%.o : %.cc
	${CC} ${CFLAGS} -c $< -o $@

%.o : %.C
	${CC} ${CFLAGS} -c $< -o $@

%.o : %.cpp
	${CC} ${CFLAGS} -c $< -o $@

%.o : %.cxx
	${CC} ${CFLAGS} -c $< -o $@

# C/C++ linker

LINKER    := g++
#LDFLAGS = 
#uncomment this when ENABLE_COMPRESSION
LDFLAGS    = -lz

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


