#!/bin/bash

bed_prefix=$1
in_dir=$2
out_dir=$3
score_files=$(ls $in_dir/*merged*)

for score_file in $score_files;
do
    filename=$(basename -- "$score_file")
    score_prefix="${filename%.*}"
    echo $score_prefix

    /wynton/group/weiss/AndrewB/plink2 --bfile $bed_prefix \
        --score $score_file 2 4 6 cols=+scoresums \
        --out "$out_dir/$score_prefix"
done
