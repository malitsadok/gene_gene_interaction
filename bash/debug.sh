#!/bin/bash
#SBATCH --time=50:00:00
#SBATCH --ntasks=2
#SBATCH --mem=100G



echo "start!!"

python /sci/labs/orzuk/mali.tsadok/UKB/code/swip.py 

echo "finish!!"

