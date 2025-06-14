HOME=.
LIB=../lib
XDRLIB=../lib/external/xdr-interface-fortran/xdrfile-1.1.4/lib
XDRINTERFACE=../lib/external/xdr-interface-fortran
EXEF=potential_of_mean_force.x
TRGT=potential_of_mean_force

# Select Fortran compiler
#
chk_ifort := $(shell which ifort    >/dev/null 2>&1 && echo 1 || echo 0)
chk_ifx   := $(shell which ifx      >/dev/null 2>&1 && echo 1 || echo 0)
chk_gfort := $(shell which gfortran >/dev/null 2>&1 && echo 1 || echo 0)


ifeq ($(FC),ifort)
	FFLAGS = -qopenmp
else ifeq ($(FC),ifx)
	FFLAGS = -qopenmp
else ifeq ($(FC),gfortran)
	FFLAGS = -fopenmp
else
	ifeq ($(chk_ifort),1)
		FC     = ifort
		FFLAGS = -qopenmp
	else ifeq ($(chk_ifx),1)
		FC     = ifx
		FFLAGS = -qopenmp
	else ifeq ($(chk_gfort),1)
		FC     = gfortran
		FFLAGS = -fopenmp
	endif
endif

# Source codes
#
OBJ=${LIB}/util.f90        \
    ${LIB}/const.f90       \
    ${LIB}/mtwist.f90      \
    ${LIB}/random.f90      \
    ${LIB}/akimautil.f90   \
    ${LIB}/grid3d.f90      \
    ${LIB}/voronoi.f90     \
    ${LIB}/input.f90       \
    ${LIB}/output.f90      \
    ${LIB}/cv.f90          \
    ${LIB}/bootstrap.f90   \
    ${HOME}/ctrl.f90       \
    ${HOME}/analyze.f90    \
    ${HOME}/main.f90

$(TRGT) : $(OBJ) 
	$(FC) $(FFLAGS) -o $(EXEF) $(OBJ) -L $(XDRLIB) -lxdrfile
	mkdir -p ../bin
	cp $(EXEF) ../bin/
	rm -f *.mod
