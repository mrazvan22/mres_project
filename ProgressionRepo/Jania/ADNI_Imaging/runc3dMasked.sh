#!/bin/bash
FILES=/Users/jaghajan/Documents/Experiments/MAPER_MASKED/Patient/*.nii
for f in $FILES
do
pref=`echo $f | sed -e 's/.nii//'`
/Users/jaghajan/Library/c3d-0.8.2-Darwin-i386/bin/c3d $f $f -label-statistics > ${pref}.txt
done


FILESC=/Users/jaghajan/Documents/Experiments/MAPER_MASKED/Control/*.nii
for f in $FILESC
do
pref=`echo $f | sed -e 's/.nii//'`
/Users/jaghajan/Library/c3d-0.8.2-Darwin-i386/bin/c3d $f $f -label-statistics > ${pref}.txt
done