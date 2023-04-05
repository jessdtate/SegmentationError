#!/bin/zsh

#./check_mapping.sh /Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Forward/vent_MD_cheat/11samples/shape_models/control_points


if [ -z $1 ]
then
echo "input directory"
exit
fi


tag=_tet2points_mapping.mat

sr_call=/Users/jess/software/SCIRun/bin/SCIRun/SCIRun_test
sr_net=/Users/jess/CIBC/FP/segmentation_error/seg_error/nets/clean_cleavermeshes.srn5
sr_net=/Users/jess/CIBC/FP/segmentation_error/seg_error/nets/clean_cleavermeshes_modes.srn5


s_files=($(ls -1 $1/model_param*.pts))

outpath=${1%/*}/meshes/
ex_path=${1%/*}/

echo $outpath
for f in ${s_files[@]}
do

filename=($(basename $f))
froot=${filename%*.pts}

#echo $f
#echo $froot

fcheck=$outpath$froot$tag

#echo $fcheck
if test -f "$fcheck";
then
    continue
else

echo $fcheck not found.  rerunning
echo $ex_path$filename
export FILENAME=$ex_path$filename
echo $sr_call -E -x -0 $sr_net
$sr_call -E -0 $sr_net
fi
done

echo "all done"
exit
