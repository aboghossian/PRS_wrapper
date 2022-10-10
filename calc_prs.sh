#!/bin/bash
# Calculates PRS scores from PRScsx output merged files
# Recommended to use SCORE1_SUM column inf follow-up (rather than AVG)

bed_prefix=$1  # bed file to score
in_dir=$2  # input directory with merged files
out_dir=$3  # output directory

# get score files from input directory
score_files=$(ls $in_dir/*merged*)

# loop through score files
for score_file in $score_files;
do
    filename=$(basename -- "$score_file")
    score_prefix="${filename%.*}"  # get file name for output
    echo $score_prefix

    # calculate scores using plink2
    /wynton/group/weiss/AndrewB/plink2 --bfile $bed_prefix \
        --score $score_file 2 4 6 cols=+scoresums \
        --out "$out_dir/$score_prefix"
done
