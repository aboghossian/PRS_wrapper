#!/bin/bash
#
#$ -S /bin/bash
#$ -o log
#$ -e log
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=20G
#$ -l h_rt=6:00:00
#$ -t 1-22

# read in arguments
while getopts ":r:b:g:p:n:d:o:" arg; do
  case $arg in
    r) # LD reference directory
      ref_dir=${OPTARG};;
    b) # bim prefix
      bim=${OPTARG};;
    g) # GWAS file (must be in PRS-CSx format)
      gwas=${OPTARG};;
    p) # population ancestry
      pop=${OPTARG};;
    n) # GWAS sample size
      n_gwas=${OPTARG};;
    d) # output directory
      out_dir=${OPTARG};;
    o) # output name (prefix)
      out_name=${OPTARG};;
  esac
done

# environment with PRS-CSx dependencies
source /wynton/home/weiss/aboghossian/PRScsx/prscsx_env/bin/activate

for phi in 1e-6 1e-4 1e-2 1;
do
  python /wynton/home/weiss/aboghossian/PRScsx/PRScsx.py \
    --ref_dir $ref_dir \
    --bim_prefix $bim \
    --sst_file $gwas \
    --pop $pop \
    --n_gwas $n_gwas \
    --out_dir $out_dir \
    --chrom $SGE_TASK_ID \
    --phi $phi \
    --out_name $out_name
done
