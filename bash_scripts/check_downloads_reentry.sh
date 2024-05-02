#!/bin/zsh

# check download folder for missing files

#./check_downloads.sh /Users/jess/CIBC/FP/segmentation_error/Dalhousie_seg/UQ_data/Forward/vent_MD_cheat/11samples/solutions/Solutions_inria_iso/inria_iso_hires_new

# TODO:  This needs to be fix since it's not working now

#drive="segmentation_uncertainty_gdrive"
#drive_path="NewCases-Solutions-hires"

drive1="SegUQ_hires_reentry_p1"
drive2="SegUQ_hires_reentry_p2"
drive_path=""

drive_cutoff=200

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

num_files_missing=0
num_files_missize=0
num_downloaded=0
download_list=()

if [ -f $tmp_filename ]
then
    echo removing file $tmp_filename
    rm $tmp_filename
fi


if [ $download -eq 1 ]
then
    printf "downloading missing files.  Takes ~40-60s a file\n"
fi


if [ $remote -eq 0 ]
then
    base="Solution_model_params_"
#    tags=( "_RVB_stim_point" "_Lvot_stim_point" "_LVV_point" "_Rvot_stim_point"  "_RVV_stim_point" )
    tags=( "" )
    ending=".mat"
    start_ind=000
#    end_ind=drive_cutoff
#    start_ind=
    end_ind=261

    (( tot_ind=end_ind+1 ))

    curr_ind=0
    #curr_file=1
    (( total_files=($end_ind-$start_ind+1)*${#tags[@]} ))

    printf "\nChecking %d files locally base on this pattern:\n" $total_files
    printf " %s%03d:%03d%s\n" $base $start_ind $end_ind $ending


    for ((k=$start_ind;k<=end_ind;k++))
    do
        need2download=0
        num=$(printf "%03d" $k)
            
        f_root=$base$num$ending
        f=$1/$f_root
        if [ $k -le $drive_cutoff ]
        then
            drive=$drive1
        else
            drive=$drive2
        fi
        remote_file=$drive:$f_root
        
#        printf "remote file: %s\n" $remote_file

        if [ ! -f $f ]
        then
            echo missing: $f_root
            need2download=1
            (( num_files_missing++ ))
        fi

        if [ $download -eq 1 ] && [ $need2download -eq 1 ]
        then
            echo downloading: $remote_file
            rclone copy --no-traverse $remote_file $1
            (( num_downloaded++ ))
        fi
    done

else

    # TODO: this isn't working.  I need to figure out how to parse this properly
    drive_contents=$(rclone ls $drive:$drive_path)
    
    IFS=$' '
    dc_array=($drive_contents)

    total_files=${#dc_array[@]}

    printf "\nChecking %d files based on the contents of the remote drive folder:\n" $total_files
    printf "%s:%s\n" $drive $drive_path

    echo "remote version not implemented, use '-l' flag"
    return
    

    if [ $check_size -eq 1 ]
    then
        printf "checking size from remote drive.\n"
    fi

    for line in ${drive_contents[@]}
    do

        rf_name=$($line | cut -f2 )

        echo $rf_name

        lf_name=$1/$rf_name

        if [ ! -f $lf_name ]
        then
            echo missing: $f_root
            need2download=1
            num_files_missing+=1
        else
            if [ $check_size -eq 1 ]
            then
                lf_size=$(wc -c <"$f")
                rf_size=$($line | cut -f1 )
                
                echo $lf_size
                echo $rf_size

                if [ $lf_size -ne $rf_size ]
                then
                    echo wrong size: $f_root
                    need2download=1
                    num_files_missize+=1
                fi
            fi
        fi

        if [ $download -eq 1 ] && [ $need2download -eq 1 ]
        then
            remote_file=$drive:$drive_path/$rf_name
            #rclone copy --no-traverse $remote_file $1
            echo $remote_file >> $tmp_filename
            num_downloaded+=1
        fi

    done

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
        printf ' --- %d mismatched file sizes found\n' $num_files_missized
    fi
fi

if [ $num_downloaded -gt 0 ]
then
    printf " --- downloaded %d files\n" $num_downloaded
    if [ $remote -eq 1 ] && [ $download -eq 1 ]
    then
      printf " --- downloading %d files\n" $num_downloaded
      rclone copy -Pv --dry-run $tmp_filename $1
      rm $tmp_filename
    fi
else
    printf " --- No files to download\n"
fi
