#!/bin/bash

# Usage: ./install.sh --compiler=<compiler type> ("gcc" or "intel")
#
# 1. If you need to install HDF5 library, please add --install-hdf5
# 2. If you wish to skip library installation, please add --skip-install-lib

# Variables
#
INSTALL_HDF5=false
SKIP_INSTALL_LIB=false
COMPILER=intel


# Parse arguments
#
while [[ $# -gt 0 ]]; do
    case "$1" in
        --install-hdf5)
            INSTALL_HDF5=true
            shift
            ;;
        --skip-install-lib)
            SKIP_INSTALL_LIB=true
            shift
            ;;
        --compiler=*)
            COMPILER="${1#*=}"
            shift
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

#
#=====================================================================
#
# Install XDR Library 
#
cwd=`pwd`
XDRPATH=$ANATRA_PATH/f90/lib/external/xdr-interface-fortran
HDFPATH=$ANATRA_PATH/f90/lib/external/hdf5 #-hdf5-1_12_3
NCPATH=$ANATRA_PATH/f90/lib/external/netcdf

if [ "$COMPILER" == "intel" ];then
  chk_ifort=`which ifort    >/dev/null 2>&1 && echo 1 || echo 0`
  chk_ifx=`which ifx        >/dev/null 2>&1 && echo 1 || echo 0`
  if [ "$chk_ifort" -eq 1 ];then
    fortcomp=ifort
  elif [ "$chk_ifx" -eq 1 ];then
    fortcomp=ifx
  fi
fi

if [ "$SKIP_INSTALL_LIB" == "true" ];then

  echo "Skip installation of libraries"

elif [ "$SKIP_INSTALL_LIB" == "false" ];then

  cd $XDRPATH
    if [ ! -e xdrfile-1.1.4 ]; then
      tar xvf xdrfile-1.1.4.tar.gz 
    fi
    cd xdrfile-1.1.4
    ./configure CC=gcc FC=gfortran --prefix=$XDRPATH/xdrfile-1.1.4
    make && make install
  cd $cwd
  
  # Install HDF5 library
  #
  if [ "$INSTALL_HDF5" == "true" ];then
    cd $HDFPATH
      if [ ! -e hdf5-hdf5-1_12_3 ];then
        tar xvf hdf5-hdf5-1_12_3.tar.gz
      fi
      cd hdf5-hdf5-1_12_3
      if [ "$COMPILER" == "gcc" ];then
        ./configure FC=gfortran CC=gcc --prefix=$HDFPATH
        make && make install
      elif [ "$COMPILER" == "intel" ];then
        ./configure FC=$fortcomp --prefix=$HDFPATH
        make && make install
      fi
    cd $cwd
  
    if [ ! -e $HDFPATH/lib/libhdf5.so ];then
      echo "Error during installation of HDF5 library to $HDFPATH/hdf5-hdf5-1_12_3/"
      exit
    fi
  
    # For building NetCDF with above HDF5
    unset PKG_CONFIG_PATH
    export CPPFLAGS="-I$HDFPATH/include"
    export LDFLAGS="-L$HDFPATH/lib"
  fi
  
  # Install NetCDF library
  #
  
  cd $NCPATH
  if [ ! -e netcdf-4.6.1 ];then
    tar xvf netcdf-4.6.1.tar.gz
  fi
  cd netcdf-4.6.1
  if [ "$COMPILER" == "gcc" ];then
    ./configure FC=gfortran CC=gcc --prefix=$NCPATH/netcdf
    make && make install
  elif [ "$COMPILER" == "intel" ];then
    ./configure FC=$fortcomp --prefix=$NCPATH/netcdf
    make && make install
  fi
  cd ..
  
  if [ ! -e $NCPATH/netcdf/lib/libnetcdf.so ];then
    echo "Error during installation of NETCDF-C library to $NCPATH/netcdf"
    exit
  fi
  
  # Install NetCDF-Fortran library
  #
  if [ ! -e netcdf-fortran-4.4.4 ];then
    tar xvf netcdf-fortran-4.4.4.tar.gz
  fi
  cd netcdf-fortran-4.4.4
  export LDFLAGS="$LDFLAGS -L$NCPATH/netcdf/lib"
  export LIBS="-lnetcdf"
  export CPPFLAGS="$CPPFLAGS -I$NCPATH/netcdf/include"
  export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$NCPATH/netcdf/lib"
  
  if [ "$COMPILER" == "gcc" ];then
    ./configure FC=gfortran CC=gcc --prefix=$NCPATH/netcdf --with-netcdf=$NCPATH/netcdf
  elif [ "$COMPILER" == "intel" ];then
    ./configure FC=$fortcomp --prefix=$NCPATH/netcdf --with-netcdf=$NCPATH/netcdf
  fi
  
  make && make install
  cd $cwd
  
  if [ ! -e $NCPATH/netcdf/lib/libnetcdff.so ];then
    echo "Error during installation of NETCDF-F library to $NCPATH/netcdf"
    exit
  fi

fi

# Install ANATRA fortran programs
#
if [ "$COMPILER" == "" ]; then
  COMPILER=intel
  echo "no compiler type is specified"
  echo ">> intel is used for compile"
elif [ "$COMPILER" == "gcc" ]||[ "$COMPILER" == "intel" ]; then
  echo "$COMPILER is used"
fi

list="center_of_mass          \
      distance                \
      lipid_order             \
      z_profile               \
      z_orient                \
      rotation                \
      radial_distr            \
      spatial_distr           \
      cluster_size            \
      extract_frame           \
      random_frame            \
      histogram               \
      moving_average          \
      average_function        \
      bootstrap_prepper       \
      iepdyn                  \
      return_prob             \
      restricted_radial       \
      potential_of_mean_force \
      interaction_energy      \
      free_volume             \
      rprate                  \
      state_define            \
      gauqm_prepper"

cwd=`pwd`
mkdir -p bin
for d in $list;do
  echo "o Installing $d ..."
  echo ""
  if [ "$COMPILER" == "gcc" ]&&[ "$d" == "en_analysis" ]; then
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

    if [ "$COMPILER" == "gcc" ]&&[ "$d" == "interaction_energy" ]; then
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

