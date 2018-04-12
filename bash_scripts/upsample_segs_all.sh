#!/bin/bash

# upsample heart segmentations

if [ -z $1 ]
then
echo need input directory
exit
fi

echo $1

files=($(ls -1 $1/*vent*nrrd | grep -v 'hires'))

for f in ${files[@]}
do

echo $f

./upsample_segs.sh $f
out=$?

done

