#!/bin/zsh

# check the subsample heart segmentations from subsample_hires.sh

#./check_hires_sampling.sh /Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Forward/vent_MD_cheat/11samples/shape_models/control_points





if [ -z $1 ]
then
echo need input directory
exit
fi



files=($(ls -1 $1/*.pts))

scirun_call=/Users/jess/software/SCIRun_testing2/bin_515/SCIRun/SCIRun_test

network=/Users/jess/CIBC/FP/segmentation_error/seg_error/nets/check_hires_subsample.srn5

for f in ${files[@]}
do

echo $f


export POINT_FILE=$f
$scirun_call -0 -E $network
out=$?

done

