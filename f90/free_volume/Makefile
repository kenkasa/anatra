FC=ifort
HOME=.
LIB=../lib
XDRLIB=../lib/external/xdr-interface-fortran/xdrfile-1.1.4/lib
XDRINTERFACE=../lib/external/xdr-interface-fortran
NCLIB=$(realpath ../lib/external/netcdf/netcdf/lib)
NCINC=$(realpath ../lib/external/netcdf/netcdf/include)
EXEF=free_volume.x
TRGT=free_volume

OBJ=${XDRINTERFACE}/xdr.F90 \
    ${LIB}/util.f90         \
    ${LIB}/const.f90        \
    ${LIB}/mtwist.f90       \
    ${LIB}/random.f90       \
    ${LIB}/grid3d.f90       \
    ${LIB}/input.f90        \
    ${LIB}/output.f90       \
    $(LIB)/netcdf.f90       \
    ${LIB}/dcd.f90          \
    ${LIB}/xtclib.f90       \
    ${LIB}/xyz.f90          \
    ${LIB}/traj.f90         \
    ${LIB}/com.f90          \
    ${LIB}/prmtop.f90       \
    ${LIB}/anaparm.f90      \
    ${LIB}/pot.f90          \
    ${HOME}/ctrl.f90        \
    ${HOME}/setup.f90       \
    ${HOME}/analyze.f90     \
    ${HOME}/main.f90

all   : ${TRGT}
clean :
	rm ${SRC}*.o *.mod

${TRGT} : ${OBJ} 
	${FC} -qopenmp -mkl -o ${EXEF} ${OBJ} -L ${XDRLIB} -lxdrfile -L $(NCLIB) -I $(NCINC) -lnetcdf -lnetcdff -Wl,-rpath,$(NCLIB)
	cp ${EXEF} ../bin/
	rm -f *.mod
