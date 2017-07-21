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


# link all .cpp .h files
#link: 
#	ln -sfr *.cpp *.h -t ${LAMMPS_SRC}
#
#unlink:
#	find ${LAMMPS_SRC} -type l -delete

# link only some files
link: 
	ln -sfr $(files) -t ${LAMMPS_SRC}

unlink:
	rm ${lammps_files}


