#!/bin/bash

# upsample heart segmentations

if [ -z $1 ]
then
echo need input directory
exit
fi

echo $1

cleaver_call=/Users/jess/software/Cleaver2/build/bin/cleaver-cli

files=($(ls -1 $1/*vent*hires*nrrd))

for f in ${files[@]}
do

echo $f
fname=${f##*/}
froot=${fname%%.*}
ft=${f%/*}

echo $cleaver_call -S -i $f -o $ft../meshes/ -n $froot -f matlab -B 1 -c 1 -j -g 0.1 -v
$cleaver_call -S -i $f -o $ft../meshes/ -n $froot -f matlab -B 1 -c 1 -j -g 0.1 -v


done

