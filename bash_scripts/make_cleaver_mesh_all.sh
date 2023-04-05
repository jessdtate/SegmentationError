#!/bin/bash

#./make_cleaver_mesh_all.sh /home/sci/jess/Segmentation_Error/Dalhousie_seg/shape_data/vent_MD/11samples/shape_models

#./make_cleaver_mesh_all.sh /Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/shape_data/vent_MD_cheat/11samples/shape_models


if [ -z $1 ]
then
echo "need to provide directory for source files"
exit
fi

tag="_*_surface_vent_vol"

s_files=($(ls -1 $1/*$tag*.nrrd))

for f in ${s_files[@]}
do

echo $f
./make_cleaver_mesh_one.sh $f
done

echo "all done"
exit
