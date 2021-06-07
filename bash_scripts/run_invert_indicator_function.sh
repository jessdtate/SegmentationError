#!/bin/bash

# ./run_invert_indicator_function.sh /Users/jess/FP/segmentation_error/Dalhousie_seg/shape_data/vent_MD_cheat/9samples/SR_surfaces


if [ -z $1 ]
then
echo need to have the directory of the source files
exit
fi

sr_call="/Users/jess/software/SCIRun/bin/SCIRun/SCIRun_test"
sr_net="/Users/jess/CIBC/FP/segmentation_error/seg_error/nets/make_inverse_indicator_funct.srn5"

tag='epi_vol'

s_files=($(ls -1 $1/*$tag*.nrrd))


for f in ${s_files[@]}
do

echo $f
export FNAME=$f
echo $sr_call -0 -x -E $sr_net
$sr_call -0 -x -E $sr_net
out=$?

echo $out

done

echo all done
exit

