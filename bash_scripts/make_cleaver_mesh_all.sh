#!/bin/bash

#./make_cleaver_mesh_all.sh 


if [ -z $1 ]
then
echo "need to provide directory for source files"
exit
fi

tag="vent_vol"

s_files=($(ls -1 $1/*$tag*.nrrd))

for f in ${s_files[@]}
do

echo $f
./make_cleaver_mesh_one.sh $f
done

echo "all done"
exit
