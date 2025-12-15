#!/bin/bash
#SBATCH --time=50:00:00
#SBATCH --ntasks=2
#SBATCH --mem=100G



echo "start !!"

python /sci/labs/orzuk/mali.tsadok/UKB/singleton_code/count_gene_from_all_gene_gene.py

echo "finish!!"

