#!/bin/zsh

# subsample heart segmentations
#./subsample_surface.sh /Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Forward/vent_MD_cheat/11samples/shape_models/control_points

if [ -z $1 ]
then
echo need input directory
exit
fi

echo $1

files=($(ls -1 $1/*.pts))

scirun_call=/Users/jess/software/SCIRun/bin_clean/SCIRun/SCIRun_test

network=/Users/jess/CIBC/FP/segmentation_error/seg_error/nets/subsample_ecgsim_solutions_newsims.srn5

for f in ${files[@]}
do

echo $f
export POINT_FILE=$f
$scirun_call -0 -E $network
out=$?

done

