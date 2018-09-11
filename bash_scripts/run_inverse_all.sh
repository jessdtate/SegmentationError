#!/bin/bash

# ./run_inverse_all.sh /Users/jess/FP/segmentation_error/Dalhousie_seg/ /Users/jess/Edgar/database/Unmodified/Human_Cardiac_Mapping/Dalhousie-2006-01-05/Interventions/BSPM/


if [ -z $2 ]
then
echo need to have two directories of the source files
exit
fi


s_files=($(ls -1 $1/*/geom/BEMatrix_python.mat))


for f in ${s_files[@]}
do

echo $f

./run_inverse_one.sh $f $2
out=$?


done

echo all done
exit

