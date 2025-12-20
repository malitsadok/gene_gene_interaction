#!/bin/bash
#!/usr/bin/env python
#SBATCH --time=50:00:00
#SBATCH --ntasks=2
#SBATCH --mem=80G




if [ "$#" -lt 1 ]; then
    echo "Please pass syn or lof"
    exit 1
fi

type=$1



subfolder_list=("${@:2}") 
csv_file="/sci/home/mali.tsadok/scripts_all_data_pipeline.csv" # Accepts additional arguments as the list of subfolders
bash_script="calc_interaction_for_all_data.sh"
	
#bash_script="calc_gene_count_service.sh"

# Extract the script, input_dir, and output_dir from the CSV file
read script input_dir output_dir < <(awk -F',' -v bash_script="$bash_script" '
NR>1 {
    gsub(/^[ \t]+|[ \t]+$/, "", $1);
    gsub(/^[ \t]+|[ \t]+$/, "", $2);
    gsub(/^[ \t]+|[ \t]+$/, "", $3);  
    gsub(/^[ \t]+|[ \t]+$/, "", $4);  
    gsub(/\r/, "", $3);  
    gsub(/\r/, "", $4);
    if ($1 == bash_script) {
        print $2 , $3, $4
        exit
    }
}' "$csv_file")

input_dir_full="${input_dir}/${type}/"
output_dir_full="${output_dir}/${type}/"
echo "Script: $script"
echo "Input Directory: $input_dir_full"
echo "Output Directory: $output_dir_full"



# If subfolder_list is empty, get all subfolders from the input_directory
input_dir_full="${input_dir}/${type}/"
output_dir_full="${output_dir}/${type}/"

if [ "$#" -eq 1 ]; then
  echo "here"
  subfolders=($(ls -d "$input_dir_full"/*/))
else
    for subfolder in "${subfolder_list[@]}"; do
      subfolders+=("$input_dir_full/$subfolder/")
    done
fi



# Loop through the subfolders in pairs
for ((i = 0; i < ${#subfolders[@]}; i++)); do
  for ((j = i + 1; j < ${#subfolders[@]}; j++)); do

    chr1=$(basename "${subfolders[$i]}")
    chr2=$(basename "${subfolders[$j]}")
	echo $chr1
	echo $chr2


	num_files1=$(find "${subfolders[$i]}" -type f | wc -l)
    num_files2=$(find "${subfolders[$j]}" -type f | wc -l)

    # Calculate the number of combinations
    number=$(($num_files1 * $num_files2))
	echo $number

	

	echo 'sbatch --array=1-$number  /sci/home/mali.tsadok/all_data/$script "$chr1" "$chr2" "$input_dir_full" "$output_dir_full"'
	
	sbatch --array=1-$number  /sci/home/mali.tsadok/all_data/$script "$chr1" "$chr2" "$input_dir_full" "$output_dir_full"

  done
done


