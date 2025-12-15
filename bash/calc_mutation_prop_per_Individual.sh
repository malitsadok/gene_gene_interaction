#!/bin/bash
#SBATCH --time=50:00:00
#SBATCH --ntasks=1              
#SBATCH --mem=80G

if [ "$#" -lt 1 ]; then
    echo "Please pass syn or lof"
    exit 1
fi

type=$1

csv_file="/sci/home/mali.tsadok/pipeline_scripts_singleton.csv"
bash_script="calc_mutation_prop_per_Individual.sh"

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

python /sci/labs/orzuk/mali.tsadok/UKB/singleton_code/$script  "$input_dir_full" "$output_dir_full"

echo "Array job submitted with $num_subfolders tasks!"
