FC=gfortran
HOME=.
LIB=../lib
EXEF=sd_gather.x
TRGT=sd_gather

OBJ=${LIB}/util.f90        \
    ${LIB}/const.f90       \
    ${LIB}/akimautil.f90   \
    ${LIB}/grid3d.f90      \
    ${LIB}/input.f90       \
    ${LIB}/output.f90      \
    ${LIB}/dcd.f90         \
    ${LIB}/traj.f90        \
    ${HOME}/sd_ctrl.f90    \
    ${HOME}/sd_gather.f90 \
    ${HOME}/sd_main.f90

all   : ${TRGT}
clean :
	rm ${SRC}*.o *.mod

${TRGT} : ${OBJ} 
	${FC} -fopenmp -o ${EXEF} ${OBJ}
	cp ${EXEF} ../bin/
