#!/bin/bash

#./cleanup_cleaver_mesh_all.sh /Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/shape_data/vent_MD_cheat/9samples/shape_models


if [ -z $1 ]
then
echo "need to provide directory for source files"
exit
fi

sr_call=/Users/jess/software/SCIRun/bin/SCIRun/SCIRun_test
sr_net=/Users/jess/CIBC/FP/segmentation_error/seg_error/nets/clean_cleavermeshes.srn5
sr_net=/Users/jess/CIBC/FP/segmentation_error/seg_error/nets/clean_cleavermeshes_modes.srn5

tag="model_param"
tag="mode_"

s_files=($(ls -1 $1/*$tag*.pts))

for f in ${s_files[@]}
do

if [[ "$f" == *"surface"* ]]; then
continue
fi

echo $f

export FILENAME=$f
echo $sr_call -E -x -0 $sr_net
$sr_call -E -x -0 $sr_net
done

echo "all done"
exit
