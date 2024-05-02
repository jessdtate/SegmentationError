#!/bin/zsh

# subsample heart segmentations
#./subsample_hires.sh /Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Forward/vent_MD_cheat/11samples/shape_models/control_points


rerun=0
expected=1

#output_dir=/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Forward/vent_MD_cheat/11samples/solutions/Solutions_inria_iso/inria_iso_hires_surface
#output_dir=/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Forward/vent_MD_cheat/11samples/solutions/Solutions_inria_iso/inria_iso_hires_LAT_sampled
output_dir=/Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Forward/vent_MD_cheat/11samples/solutions/Solutions_inria_iso/inria_iso_reentry_hires_sampled


while getopts "c" flag
do
    case "${flag}" in
        c)
        rerun=1
        ;;
        \?)
        echo "Invalid option: -$OPTARG" >&2
        exit 1
        ;;
    esac
done

shift $((OPTIND - 1))

echo $rerun

if [ -z $1 ]
then
echo need input directory
exit
fi


echo $1
if [ $rerun -eq 0 ]
then
echo $output_dir
fi

files=($(ls -1 $1/*.pts))

scirun_call=/Users/jess/software/SCIRun_testing2/bin_headless/SCIRun/SCIRun_test
#scirun_call=/Users/jess/software/SCIRun_testing2/bin_515/SCIRun/SCIRun_test

#network=/Users/jess/CIBC/FP/segmentation_error/seg_error/nets/subsample_hires_solutions.srn5
network=/Users/jess/CIBC/FP/segmentation_error/seg_error/nets/subsample_hires_solutions_reentry.srn5
#network=/Users/jess/CIBC/FP/segmentation_error/seg_error/nets/subsample_hires_LAT.srn5

for f in ${files[@]}
do

echo $f


if [ $rerun -eq 0 ]
then
f_path="$(dirname "$f")"
f_base="$(basename "$f" .pts)"
#out_f=$output_dir"/Solution_"$f_base*surface.mat
#ofs=($(ls -1 $output_dir"/Solution_"$f_base*"surface.mat"))

out_f=$output_dir"/Solution_"$f_base*_sampled.mat
ofs=($(ls -1 $output_dir"/Solution_"$f_base*"_sampled.mat"))

#echo ${#ofs[@]}
#echo $expected
if [ ${#ofs[@]} -ge $expected ]
then
echo "${#ofs[@]}"
continue
#else
#echo "running"
fi

fi


export POINT_FILE=$f
$scirun_call -0 -E $network
out=$?

done

