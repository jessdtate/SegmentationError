#!/bin/zsh

# subsample heart segmentations
#./fix_matfile.sh /Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Forward/vent_MD_cheat/11samples/solutions/Solutions_ECGsim/activation_times_sampled

if [ -z $1 ]
then
echo need input directory
exit
fi

echo $1

files=($(ls -1 $1))

scirun_call=/Users/jess/software/SCIRun_testing2/bin_headless/SCIRun/SCIRun_test

network=/Users/jess/CIBC/FP/segmentation_error/seg_error/nets/fix_matfile.srn5

for f in ${files[@]}
do

echo $f
export FILENAME=$1/$f
$scirun_call -0 -x -E $network
out=$?

done

