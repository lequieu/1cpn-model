# ------------------------------------------------------------
# soft link  1CPN src files to LAMMPS_SRC directory
# ------------------------------------------------------------

# can set LAMMPS_SRC from command line using "make LAMMPS_SRC=MYPATH"

LAMMPS_SRC ?= '../../../lammps-1cpn/src/'

files = angle_orient.cpp \
		angle_orient.h \
		angle_orient_cosine.cpp \
		angle_orient_cosine.h \
        angle_wlc_twist.cpp \
        angle_wlc_twist.h \
        pair_gauss_aniso.cpp \
        pair_gauss_aniso.h \
        fix_restrain_orient.cpp \
        fix_restrain_orient.h \
        pair_zewdie.cpp\
        pair_zewdie.h

lammps_files := $(addprefix $(LAMMPS_SRC)/,$(files))
local_files := $(addprefix $(PWD)/,$(files))


# link all .cpp .h files
#link: 
#	ln -sfr *.cpp *.h -t ${LAMMPS_SRC}
#
#unlink:
#	find ${LAMMPS_SRC} -type l -delete


help:
	@echo "Please use one of the following options:"
	@echo "  make link       to soft link the 1CPN-model C++ files into the LAMMPS src directory "
	@echo "  make unlink     to remove the soft links between the 1CPN-model and LAMMPS "
	@echo ""
	@echo "Note: Path to LAMMPS source is  LAMMPS_SRC=${LAMMPS_SRC}"
	@echo "This can be changed using 'make link LAMMPS_SRC=<new path>.'"

# link only some files
link: 
	#ln -sfr $(files) -t ${LAMMPS_SRC} # Older ln doesnt permit the -r flag
	ln -sf $(local_files) -t ${LAMMPS_SRC}

unlink:
	rm ${lammps_files}


