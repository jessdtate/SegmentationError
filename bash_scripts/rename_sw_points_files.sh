#!/bin/bash

# upsample heart segmentations

if [ -z $1 ]
then
data_dir=$PWD
$tag=".isores.center.pad.cropped.tpSmoothDT_world"
echo using current directory and tag $tag
elif [ -z $2 ]
then
echo using current directory
data_dir=$PWD
tag=$1

else
data_dir=$1
tag=$2
fi


echo using directory $data_dir
echo using tag $tag


files=($(ls -1 $data_dir/*$tag*.particles))

for f in ${files[@]}
do


#fname=${f##*/}
#froot=${fname%%.*}
#ft=${f%/*}

f_out=$f".txt"

cp $f $f_out

echo $f
echo "renamed to "
echo $f_out

done

