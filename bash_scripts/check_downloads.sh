#!/bin/zsh

# check download folder for missing files

#./check_downloads.sh /Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Forward/vent_MD_cheat/11samples/solutions/Solutions_inria_iso/inria_iso_hires_new

# setting up a drive: `rclone config`
# client ID found on google developer console, API & services, credentials
# https://console.cloud.google.com/apis/credentials?project=segmentationuq
# requires authentification to add a new drive
#
# list available drives with:
# `rclone listremotes`

drive="segmentation_uncertainty_gdrive"
new_drive=1

python_exe=python3.10

# segmentation_uncertainty_gdrive:SimulationResults/NZ_Eikonal/NewCases-Solutions-hires
if [ $new_drive -eq 1 ]
then
  drive_path="SimulationResults/NZ_Eikonal/NewCases-Solutions-hires"
  #local only
  tags=( "_RVB_stim_point" "_Lvot_stim_point" "_LVV_point" "_Rvot_stim_point"  "_RVV_stim_point" )
else
  drive_path="SimulationResults/NZ_Eikonal/solutions_hires"
  #local only
  tags=( "_sinus_stim_points" "_LV_stim_point" "_RV_stim_point" "_septal_stim_point" "_apex_stim_point" )
fi

# local only.  file expectations
base="Solution_model_params_"
ending=".mat"
start_ind=000
end_ind=261



download=0
check_size=0
remote=1

json_script=../py_scripts/jsonfilesize.py

tmp_filename=tmp_download_list.txt

while getopts "lsd" flag
do
    case "${flag}" in
        d)
            download=1
        ;;
        l)
            remote=0
        ;;
        s)
            check_size=1
            remote=1
        ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
        ;;
    esac
done

shift $((OPTIND - 1))


if [ -z $1 ]
then
    echo need input directory
    exit
fi

need2download=0

let num_files_missing=0
let num_files_missize=0
let num_downloaded=0
download_list=()

if [ -f $tmp_filename ]
then
    echo removing file $tmp_filename
    rm $tmp_filename
fi


if [ $download -eq 1 ]
then
    printf "downloading missing files.  Takes ~40-60s a file on a good network \n"
fi


if [ $remote -eq 0 ]
then
    (( tot_ind=end_ind+1 ))

    curr_ind=0
    #curr_file=1
    (( total_files=($end_ind-$start_ind+1)*${#tags[@]} ))

    printf "\nChecking %d files locally base on this pattern:\n" $total_files
    printf " %s%03d:%03d_[tag]%s\n" $base $start_ind $end_ind $ending


    for ((k=$start_ind;k<=end_ind;k++))
    do
        num=$(printf "%03d" $k)
        for t in ${tags[@]}
        do
            need2download=0
            f_root=$base$num$t$ending
            f=$1/$f_root

            remote_file=$drive:$drive_path/$f_root

            if [ ! -f $f ]
            then
                echo missing: $f_root
                need2download=1
                num_files_missing+=1
            fi

            if [ $download -eq 1 ] && [ $need2download -eq 1 ]
            then
                echo downloading: $f_root
                rclone copy --no-traverse $remote_file $1
                num_downloaded+=1
            fi

        done
    done

else

    printf "Querying remote folder: %s:%s\n" $drive $drive_path
    
#    drive_contents=$(rclone ls $drive:$drive_path)
#    echo $drive_contents
#    IFS=$' '
#    dc_array=($drive_contents)
#
#    total_files=${#dc_array[@]}
    
    

    size_json=$(rclone size --json $drive:$drive_path)
    
    num_files=$($python_exe $json_script $size_json count)
    
#    echo $num_files

    printf "\nChecking %d files based on the contents of the remote drive folder:\n" $num_files
    

#    echo "remote version not implemented, use '-l' flag"
#    return
    
    if [ $check_size -eq 1 ]
    then
        printf "checking size from remote drive.\n"
    fi
    
    while read -r size name; do
#        echo $size
#        echo $name
        need2download=0
        lf_name=$1/$name
        if [ ! -f $lf_name ]
        then
            echo missing: $name
            need2download=1
            num_files_missing+=1
        else
            if [ $check_size -eq 1 ]
            then
                lf_size=$(wc -c <"$lf_name")
                
                

                if [ $lf_size -ne $size ]
                then
                    echo wrong size: $lf_name
                    echo size $lf_size vs $size
                    need2download=1
                    num_files_missize+=1
                fi
            fi
        fi
        
        if [ $download -eq 1 ] && [ $need2download -eq 1 ]
        then
#            remote_file=$drive:$drive_path/$rf_name
            echo dowloading: $name $size
            
            rclone copy --no-traverse $drive:$drive_path/$name $1
            echo $name >> $tmp_filename
            num_downloaded+=1
        fi
        
    done < <(rclone ls $drive:$drive_path)


#    for line in ${total_files[@]}
#    do
##        echo $line
#
#        rf_name=$($line | cut -f2 )
#
##        echo $rf_name
#
#        lf_name=$1/$rf_name
#
#        if [ ! -f $lf_name ]
#        then
#            echo missing: $rf_name
#            need2download=1
#            num_files_missing+=1
#        else
#            if [ $check_size -eq 1 ]
#            then
#                lf_size=$(wc -c <"$f")
#                rf_size=$($line | cut -f1 )
#
#                echo $lf_size
#                echo $rf_size
#
#                if [ $lf_size -ne $rf_size ]
#                then
#                    echo wrong size: $f_root
#                    need2download=1
#                    num_files_missize+=1
#                fi
#            fi
#        fi
#
#        if [ $download -eq 1 ] && [ $need2download -eq 1 ]
#        then
#            remote_file=$drive:$drive_path/$rf_name
#            #rclone copy --no-traverse $remote_file $1
#            echo $remote_file >> $tmp_filename
#            num_downloaded+=1
#        fi
#
#    done

fi


printf '\nFinished!\n'

printf '\nReport:\n'
printf '================\n'

if [ $num_files_missing -eq 0 ]
then
    printf ' --- No missing files\n'
else
    printf ' --- %d missing files\n' $num_files_missing
fi


if [ $check_size -eq 1 ]
then
    if [ $num_files_missize -eq 0 ]
    then
        printf ' --- No mismatched file sizes found\n'
    else
        printf ' --- %d mismatched file sizes found\n' $num_files_missize
    fi
fi

if [ $num_downloaded -gt 0 ]
then
    printf " --- download %d files\n" $num_downloaded
#    rclone copy -Pv --from-file $tmp_filename $1
    rm $tmp_filename
elif [ $download -eq 1 ]
then
    printf " --- No files to download\n"
fi
