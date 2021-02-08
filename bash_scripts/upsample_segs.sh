#!/bin/bash

# upsample heart segmentations

if [ -z $1 ]
then
echo need 1 imput input filename
exit
fi

fname=${1##*/}
froot=${fname%%.*}
ft=${1%/*}

#echo variables:
#echo fname = $fname
#echo froot = $froot
#echo ft = $ft

seg3d_bin=/Users/jess/software/Seg3D/bin/Seg3D/
seg3d_call=Seg3D2.app/Contents/MacOS/Seg3D2
curr_dir=$PWD
s_file=$ft'/'$froot'_tmp.py'

echo -e "exec(open('/Users/jess/FP/segmentation_error/seg_error/py_scripts/upsample_segs.py').read()) \nupsample_heart('"$1"', '"$ft'/'"','"$froot'_hires'"')" > $s_file

echo $seg3d_call --python=$s_file
cd $seg3d_bin
#$seg3d_call --headless --python=$s_file
$seg3d_call --python=$s_file
out=$?
cd $curr_dir

rm $s_file
