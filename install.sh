#!/bin/bash

# Usage: ./install.sh --compiler=<compiler type> ("gcc" or "intel")
#
# 1. If you need to install HDF5 library, please add --install-hdf5
# 2. If you wish to skip library installation, please add --skip-install-lib

# Variables
#
INSTALL_HDF5=""
SKIP_INSTALL_LIB=""
COMPILER=intel

# Parse arguments
#
while [[ $# -gt 0 ]]; do
    case "$1" in
        --install-hdf5)
            INSTALL_HDF5="--install-hdf5"
            shift
            ;;
        --skip-install-lib)
            SKIP_INSTALL_LIB="--skip-install-lib"
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
if [ "$COMPILER" == "" ]; then
  COMPILER=intel
  echo "no compiler type is specified"
  echo ">> intel is used for compile"
elif [ "$COMPILER" == "gcc" ]||[ "$COMPILER" == "intel" ]; then
  echo "$COMPILER is used"
fi

echo "#############################################"
echo "      Installation of ANATRA"
echo "#############################################"
echo ""


echo "---------------------------------------------"
echo "Step 1. Setup ANATRA_PATH variable"
echo "---------------------------------------------"
echo ""

cwd=`pwd`
var="export ANATRA_PATH=$cwd"
var2="export PATH=\$PATH:\$ANATRA_PATH/bin"
var3="export LD_LIBRARY_PATH=\$ANATRA_PATH/f90/lib/external/netcdf/netcdf/lib:\$LD_LIBRARY_PATH"

use_bash=false
use_zsh=false
if   [ "$(basename "$SHELL")" == "bash" ]; then
  use_bash=true
elif [ "$(basename "$SHELL")" == "zsh" ]; then
  use_zsh=true
fi

if   [ "$use_bash" == "true" ]; then
  chk=`grep "ANATRA_PATH" ~/.bashrc | wc -l`
  if [ $chk -eq 0 ];then
    echo "$var"  >> ~/.bashrc
    echo "$var2" >> ~/.bashrc
    echo "$var3" >> ~/.bashrc
    echo "ANATRA_PATH has been defined in ~/.bashrc"
  else
    echo "ANATRA_PATH is already set in ~/.bashrc"
  fi 
elif [ "$use_zsh"  == "true" ]; then 
  chk=`grep "ANATRA_PATH" ~/.zshrc | wc -l`
  if [ $chk -eq 0 ];then
    echo "$var"  >> ~/.zshrc
    echo "$var2" >> ~/.zshrc
    echo "$var3" >> ~/.zshrc
    echo "ANATRA_PATH has been defined in ~/.zshrc"
  else
    echo "ANATRA_PATH is already set in ~/.zshrc"
  fi 
fi

if [ "$use_bash" == "true" ];then
  source ~/.bashrc
elif [ "$use_zsh"  == "true" ]&&[ $chnk -eq 0 ]; then
  echo "Please rerun ./install.sh after running 'source ~/.zshrc'"
  exit
fi	

echo ""
echo "---------------------------------------------"
echo "Step 2. Compiling ANATRA Fortran programs"
echo "---------------------------------------------"
echo ""

cwd=`pwd`

cd f90 
./install.sh --compiler=$COMPILER $INSTALL_HDF5 $SKIP_INSTALL_LIB
cd $cwd

cd utility
./install.sh --compiler=$COMPILER
cd $cwd


