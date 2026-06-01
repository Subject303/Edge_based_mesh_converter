#!/bin/bash
set -e

rm -f "../modfiles/*"
rm -f "../converter"

LOG="./logfile.txt"
rm -f $LOG
touch $LOG
chmod 777 $LOG

COMPLIST=" ../src/src_shared/* "

COMPLIST=$COMPLIST" ../src/src_serial/*"

COMP_OPTIONS=" -fdec-format-defaults -ffree-line-length-none -ffpe-trap=zero,overflow  "

gfortran $COMP_OPTIONS -o "../converter" $COMPLIST -J"../modfiles" | $LOG

if [ -f "../converter"  ]; then
    echo " passed compilation "
else
    echo "compilation failed"
    exit
fi

echo ' running python preprocessor '

python3 "../python_scripts/case_preprocessor.py"

echo ' running main converter program serial '

../converter | $LOG

echo ' finished. '
