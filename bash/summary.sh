#!/bin/bash
#SBATCH --time=50:00:00
#SBATCH --ntasks=2
#SBATCH --mem=256G



echo "start!!"

python /sci/labs/orzuk/mali.tsadok/UKB/singleton_code/summary_for_all_SNV.py

echo "finish!!"

