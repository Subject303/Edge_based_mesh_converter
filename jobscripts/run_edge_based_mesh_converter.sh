#!/bin/bash
set -e

rm -f "../modfiles/*"
rm -f "../converter"

LOG="logfile.txt"

COMPLIST=" ../src/src_shared/* "

COMPLIST=$COMPLIST" ../src/src_serial/*"

COMP_OPTIONS=" -fdec-format-defaults -ffree-line-length-none -ffpe-trap=zero,overflow  "

gfortran $COMP_OPTIONS -o "../converter" $COMPLIST -J"../modfiles" &> $LOG

if [ -f "../converter"  ]; then
    echo " passed compilation " &>> $LOG
else
    echo "compilation failed" &>> $LOG
    exit
fi

echo ' running python preprocessor ' &>> $LOG

python3 "../python_scripts/case_preprocessor.py" &>> $LOG

echo ' running main converter program serial ' &>> $LOG

../converter &>> $LOG

echo ' finished. ' &>> $LOG
