#!/bin/bash


if [ -z $2 ]
then
echo need to input BEMatrix and BSPM folder
exit
fi



scirun_call=/Users/jess/software/SCIRun/bin/SCIRun/SCIRun_test
#s_file=/Users/jess/FP/Forward_Error/forward_error/py_files/run_evaluate.py
scirun_net=/Users/jess/FP/segmentation_error/seg_error/nets/inverse_solution.srn5


echo $1
echo $2

fname=${1##*/}
froot=${fname%%.*}
ft=${1%/*}

s_file=$ft'/'$froot'_tmp.py'

echo -e "import threading\nimport os\nimport fnmatch\nimport sys\nsource_file='"$1"'\ndirname = '"$2"'\nfiles = os.listdir(dirname)\nscirun_load_network('"$scirun_net"')\nscirun_set_module_state('ReadMatrix:1','Filename',source_file) \nfor f in files:\n  print(f)\n  scirun_set_module_state('ReadMatrix:0','Filename',dirname + f)\n  scirun_execute_all()\nscirun_force_quit()" > $s_file


echo running $1
echo $scirun_call -0 -S $s_file
$scirun_call -0 -x -S $s_file
out=$?
echo $out


rm $s_file

exit 0
