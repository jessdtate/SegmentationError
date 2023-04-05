#!/bin/bash

# mesh_segmentations with cleaver

if [ -z $1 ]
then
echo need input directory
exit
fi

echo $1

cleaver_bin=/home/sci/jess/cibcsm_builds/Cleaver2/build/bin
cleaver_bin=/Users/jess/software/Cleaver2/build/bin
cleaver_call=$cleaver_bin/cleaver-cli

tag=ventricles_hires

files=($(ls -1 $1/*$tag*nrrd))

for f in ${files[@]}
do

echo $f

in_dir=$1
filename=($(basename $f))

froot=${filename%$tag*}
echo $froot

out_dir=$in_dir/cleaver_meshes/

F=1
B=.7

outname=$filename"_MD_F"$F

echo $cleaver_call -I -i $in_dir/$filename  -j -e -v -f matlab -o $out_dir -B $B -F $F -n $outname

$cleaver_call -I -i $in_dir/$filename -j -e -v -f matlab -o $out_dir -B $B -F $F -n $outname

done

