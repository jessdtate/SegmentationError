#!/bin/bash

# upsample heart segmentations

if [ -z $1 ]
then
echo need 1 imput input filename
exit
fi

scirun_call=/Users/jess/software/SCIRun/bin/SCIRun/SCIRun_test
scirun_net=/Users/jess/FP/segmentation_error/seg_error/nets/make_forward_matrix.srn5


echo $1

files=($(ls -1 $1/*/geom/heart_surface_small.stl))

for f in ${files[@]}
do

echo $f

export HEART_SURFACE=$f

echo $scirun_call -E $scirun_net
$scirun_call -E -x -0 $scirun_net
out=$?

done
