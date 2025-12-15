#!/bin/bash
#!/usr/bin/env python
#SBATCH --time=50:00:00
#SBATCH --ntasks=2
#SBATCH --mem=200G

chr1=$1
chr2=$2
path=$3
result_dir=$4

path1="$path$chr1"
echo $path1
path2="$path$chr2"

echo $path2
main_dir="$path"

echo $result_dir

# Populate list1 with all the file names from path1 without the .csv extension
list1=()
while IFS= read -r -d '' file; do
  filename=$(basename "$file" .csv)  # Extract filename without the extension
  list1+=("$filename")
done < <(find "$path1" -type f -name "*.csv" -print0)

# Populate list2 with all the file names from path2 without the .csv extension
list2=()
while IFS= read -r -d '' file; do
  filename=$(basename "$file" .csv)  # Extract file name without the path
  list2+=("$filename")
done < <(find "$path2" -type f -name "*.csv" -print0)

# Calculate the total number of combinations
total_combinations=$(( ${#list1[@]} * ${#list2[@]} ))

# Calculate the index of the current combination
index=$(( ($SLURM_ARRAY_TASK_ID - 1) % $total_combinations ))

# Get the parameters for this task
param1_index=$(( $index / ${#list2[@]} ))
param2_index=$(( $index % ${#list2[@]} ))
param1=${list1[$param1_index]}
param2=${list2[$param2_index]}

echo $param1
echo $param2



python /sci/labs/orzuk/mali.tsadok/UKB/singleton_code/genes_interactions_service.py "$param1" "$param2"  "$path" "$result_dir"


