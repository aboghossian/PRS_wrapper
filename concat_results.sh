#!/bin/bash
# Concatenates PRS-CSx results by phi value (single ancestry)
in_dir=$1  # input directory
phis=$2  # phi values (separated by comma)

phis_split=$(echo $phis | tr "," "\n")  # split into individual values

# loop through phis
for phi in $phis_split; do
    echo $phi

    # concatenate files with that value in name
    cat $in_dir/*pst\_eff*$phi*.txt > $in_dir/$phi\_merged.txt
done
