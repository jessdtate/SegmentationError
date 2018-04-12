#!/bin/bash

# upsample heart segmentations

if [ -z $1 ]
then
echo need input directory
exit
fi

echo $1

files=($(ls -1 $1/*torso*.nrrd | grep -v 'points'))

for f in ${files[@]}
do

echo $f

./torso_points_one.sh $f
out=$?

done

