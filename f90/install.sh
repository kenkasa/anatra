#!/bin/bash

# Usage: ./install.sh <compiler type> ("gcc" or "intel")
#
compiler=$1   # "intel" or "gcc"
#
#=====================================================================
#
# ... install xdr library
#
cwd=`pwd`
XDRPATH=$ANATRA_PATH/f90/lib/external/xdr-interface-fortran
NCPATH=$ANATRA_PATH/f90/lib/external/netcdf

if [ "$compiler" != "fugaku" ]; then
  cd $XDRPATH
    if [ ! -e xdrfile-1.1.4 ]; then
      tar xvf xdrfile-1.1.4.tar.gz 
    fi
    cd xdrfile-1.1.4
    ./configure CC=gcc FC=gfortran --prefix=$XDRPATH/xdrfile-1.1.4
    make && make install
  cd $cwd

  cd $NCPATH
  if [ ! -e netcdf-4.6.1 ];then
    tar xvf netcdf-4.6.1.tar.gz
  fi
  cd netcdf-4.6.1
  if [ "$compiler" == "gcc" ];then
    ./configure FC=gfortran CC=gcc --prefix=$NCPATH/netcdf
    make && make install
  elif [ "$compiler" == "intel" ];then
    ./configure FC=ifort CC=icc --prefix=$NCPATH/netcdf
    make && make install
  fi
  cd ..

  if [ ! -e netcdf-fortran-4.4.4 ];then
    tar xvf netcdf-fortran-4.4.4.tar.gz
  fi
  cd netcdf-fortran-4.4.4
  export LDFLAGS="-L$NCPATH/netcdf/lib"
  export LIBS="-lnetcdf"
  export CPPFLAGS="-I$NCPATH/netcdf/include"
  export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$NCPATH/netcdf/lib"

  if [ "$compiler" == "gcc" ];then
    ./configure FC=gfortran CC=gcc --prefix=$NCPATH/netcdf --with-netcdf=$NCPATH/netcdf
  elif [ "$compiler" == "intel" ];then
    ./configure FC=ifort CC=icc --prefix=$NCPATH/netcdf --with-netcdf=$NCPATH/netcdf
  fi

  make && make install
  cd $cwd
fi
#
# ... install ANATRA fortran programs
#
if [ "$compiler" == "" ]; then
  compiler=intel
  echo "no compiler type is specified"
  echo ">> intel is used for compile"
elif [ "$compiler" == "gcc" ]||[ "$compiler" == "intel" ]; then
  echo "$compiler is used"
fi

list="center_of_mass          \
      lipid_order             \
      z_profile               \
      z_orient                \
      radial_distr            \
      spatial_distr           \
      extract_frame           \
      random_frame            \
      histogram               \
      average_function        \
      trans_prob              \
      return_prob             \
      restricted_radial       \
      potential_of_mean_force \
      interaction_energy      \
      state_define"

cwd=`pwd`
mkdir -p bin
for d in $list;do
  echo "o Installing $d ..."
  echo ""
  if [ "$compiler" == "gcc" ]&&[ "$d" == "en_analysis" ]; then
    echo "Compiler: gcc  Program: en_analysis"
    echo ">> Skipped"
    echo ""
    continue
  fi
  cd $d

  make -f Makefile

  cd $cwd 
  echo ">> Finished"
  echo "" 
done

chk=0
for d in $list;do
  if [ ! -e ./bin/${d}.x ];then

    if [ "$compiler" == "gcc" ]&&[ "$d" == "interaction_energy" ]; then
      continue
    fi

    echo "Installation of $d is failed."
    echo "Please contact the developers"
    echo "if the problem is due to bugs."
    echo ""
    chk=`expr $chk + 1` 
  fi 
done

if [ $chk -eq 0 ];then
  echo "-------------------------------------------------"
  echo "Installation of ANATRA fortran programs have been"
  echo "succesfully finished!!"
  echo "-------------------------------------------------"
else
  echo "-------------------------------------------------"
  echo "$chk errors occured during the installation."
  echo "Installation terminated abnormally."
  echo "-------------------------------------------------"
fi

