#!/bin/bash

# ./make_seg3d_ind_function_MD.sh /Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/shape_data/vent_MD_cheat/9samples/shape_models/model_params_001.pts

if [ -z $1 ]
then
echo need point file to run
exit
fi

fname=${1##*/}
froot=${fname%.*}
ft=${1%/*}

s_file=$ft'/'$froot'_tmp.py'

seg3d_call='/Applications/Seg3D2.app/Contents/MacOS/Seg3D2'
script_file='/Users/jess/CIBC/FP/segmentation_error/seg_error/py_scripts/implicit_function_From_points_multidomain.py'

vol_path="/Users/jess/software/ShapeWorks/Examples/Python/Test_epi_centered/PrepOutput/groom_and_meshes/"
vol_name='BO_epi_centered.isores.center.pad.cropped.ISO.nrrd'

vol_fname=$vol_path$vol_name
out_fname=$froot"_surface"



echo -e "import os\nsys.path.append('/Library/Frameworks/Python.framework/Versions/3.4/lib/python3.4/site-packages/')\nexec(open('"$script_file"').read())\npoints = load_points('"$1"')\nimplicit_model('"$vol_fname"', points , '"$ft/"', '"$out_fname"')" > $s_file

echo running $1
echo $seg3d_call --headless --python=$s_file
$seg3d_call --headless --python=$s_file
out=$?
echo $out

rm $s_file

exit 0

