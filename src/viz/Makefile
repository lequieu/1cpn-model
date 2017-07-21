CC =     g++
#CFLGS      =  -O2 -finline-functions -funroll-loops
LAMMPS_SRC ?= '../../../lammps-1cpn/src/'
BINDIR = '../../bin'
CFLGS      =  -ggdb 
LIBRY      =  -L/usr/local/lib -lm 
INCLD      = -I../include -I${LAMMPS_SRC}

all: dump_to_xyz dump_to_xyz_fvu

dump_to_xyz : dump_to_xyz.cpp 
	${CC} ${CFLGS} ${INCLD} dump_to_xyz.cpp -o  ${BINDIR}/dump-to-xyz

dump_to_xyz_fvu : dump_to_xyz_fvu.cpp 
	${CC} ${CFLGS} ${INCLD} dump_to_xyz_fvu.cpp -o  ${BINDIR}/dump-to-xyz-fvu
