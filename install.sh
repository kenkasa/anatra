#!/bin/bash

# Usage: ./install.sh <compiler type> ("gcc" or "intel")
#
compiler=$1   # "intel" or "gcc"

#
#=====================================================================
#
if [ "$compiler" == "" ]; then
  compiler=intel
  echo "no compiler type is specified"
  echo ">> gcc is used for compile"
elif [ "$compiler" == "gcc" ]||[ "$compiler" == "intel" ]||[ "$compiler" == "fugaku" ]; then
  echo "$compiler is used"
fi

list="tcl/trj_convert/tr_convert.tcl        \
      tcl/comcrd_analysis/cm_analysis.tcl   \
      tcl/bangle_analysis/ba_analysis.tcl   \
      tcl/distance_analysis/ds_analysis.tcl \
      tcl/mindist_analysis/mi_analysis.tcl  \
      tcl/scd_analysis/sc_analysis.tcl      \
      tcl/zprof_analysis/zp_analysis.tcl    \
      tcl/rdf_analysis/rd_analysis.tcl      \
      tcl/dipole_analysis/dp_analysis.tcl   \
      tcl/pbcdist_analysis/pd_analysis.tcl"


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
chk=`grep "ANATRA_PATH" ~/.bashrc | wc -l`

if [ $chk -eq 0 ]; then
  echo "$var"  >> ~/.bashrc
  echo "$var2" >> ~/.bashrc
  echo "$var3" >> ~/.bashrc
  echo "ANATRA_PATH has been defined in ~/.bashrc as"
  echo "$ANATRA_PATH"
else
  echo "ANATRA_PATH is already defined in ~/.bashrc as"
  echo "$ANATRA_PATH"
fi

source ~/.bashrc

echo ""
echo "---------------------------------------------"
echo "Step 2. Compiling ANATRA Fortran programs"
echo "---------------------------------------------"
echo ""

cwd=`pwd`

cd f90 
./install.sh $compiler
cd $cwd

#if [ "$compiler" != "fugaku" ];then
#  cd utility
#  ./install.sh $compiler
#  cd $cwd
#fi
