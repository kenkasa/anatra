#!/bin/bash

# Usage: ./install.sh --compiler=<compiler type> ("gcc" or "intel")
#

# Variables
#
COMPILER=intel

# Parse arguments
#
while [[ $# -gt 0 ]]; do
    case "$1" in
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


#=====================================================================
#
# ... install ANATRA fortran programs
#
if [ "$COMPILER" == "" ]; then
  COMPILER=intel
  echo "no compiler type is specified"
  echo ">> intel is used for compile"
elif [ "$COMPILER" == "gcc" ]||[ "$COMPILER" == "intel" ]; then
  echo "$COMPILER is used"
fi

list="catcrd"

cwd=`pwd`
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

