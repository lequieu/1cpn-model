CC =     g++
LAMMPS_SRC ?= '../../../lammps-1cpn/src/'
BINDIR = '../../bin'
CFLGS      =  -ggdb -std=c++11
LIBRY      =  -L/usr/local/lib -lm 
INCLD      = -I../include -I${LAMMPS_SRC}

DEBUG ?= 0
ifeq ($(DEBUG),1)
	CFLGS =  -ggdb -Wall -std=c++11
else
	CFLGS =  -O2 -finline-functions -funroll-loops -std=c++11
endif

all: dump_to_xyz 

dump_to_xyz : dump_to_xyz.cpp 
	${CC} ${CFLGS} ${INCLD} dump_to_xyz.cpp -o  ${BINDIR}/dump-to-xyz

#dump_to_xyz_fvu : dump_to_xyz_fvu.cpp 
#	${CC} ${CFLGS} ${INCLD} dump_to_xyz_fvu.cpp -o  ${BINDIR}/dump-to-xyz-fvu
