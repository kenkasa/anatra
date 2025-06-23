#!/bin/bash

# Usage: ./install.sh <compiler type> ("gcc" or "intel")
#
compiler=$1   # "intel" or "gcc"
#
#=====================================================================
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

list="catcrd"

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
  if [ ! -e ../bin/${d} ];then

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

