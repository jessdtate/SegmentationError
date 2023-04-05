#!/bin/bash

#./make_cleaver_mesh_one.sh model_params_261_surface_vent_vol.nrrd


if [ -z $1 ]
then
echo "need to provide 'vent_vol' file"
exit
fi

if [[ "$1" != *"vent_vol"* ]]; then
  echo "filename is the wrong format provide 'vent_vol' file"
  exit
fi


cleaver_bin=/home/sci/jess/cibcsm_builds/Cleaver2/build/bin
cleaver_bin=/Users/jess/software/Cleaver2/build/bin
cleaver_call=$cleaver_bin/cleaver-cli

#out_dir=/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/shape_data/vent_MD_cheat/9samples/shape_models/cleaver_meshes/

in_dir=($(dirname $1))
filename=($(basename $1))
froot=${filename%_vent_vol*}
echo $froot

out_dir=$in_dir/cleaver_meshes/

lv_file=$froot"_lv_vol.nrrd"
rv_file=$froot"_rv_vol.nrrd"
air_file=$froot"_air_vol.nrrd"


F=1
B=.7

outname=$filename"_MD_F"$F

echo $cleaver_call -I -i $in_dir/$filename  -j -e -v -f matlab -o $out_dir -B $B -F $F -n $outname

$cleaver_call -I -i $in_dir/$filename -j -e -v -f matlab -o $out_dir -B $B -F $F -n $outname

#echo $cleaver_call -I -i $in_dir/$filename $in_dir/$lv_file $in_dir/$rv_file $in_dir/$air_file -j -e -v -f matlab -o $out_dir -B $B -F $F -n $outname
#
#$cleaver_call -I -i $in_dir/$filename $in_dir/$lv_file $in_dir/$rv_file $in_dir/$air_file -j -e -v -f matlab -o $out_dir -B $B -F $F -n $outname
