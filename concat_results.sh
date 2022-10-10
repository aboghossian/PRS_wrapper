#!/bin/bash
in_dir=$1  # input directory
phis=$2  # phi values (separated by comma)

phis_split=$(echo $phis | tr "," "\n")
for phi in $phis_split; do
    echo $phi
    cat "$in_dir/*pst\_eff*$phi*.txt" > "$in_dir/$phi\_merged.txt"
done
