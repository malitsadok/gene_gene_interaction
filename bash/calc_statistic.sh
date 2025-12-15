#!/bin/bash
#SBATCH --time=50:00:00
#SBATCH --ntasks=1              
#SBATCH --mem=150G

if [ "$#" -lt 1 ]; then
    echo "Please pass syn or lof"
    exit 1
fi

type=$1
subfolder_list=("${@:2}")  # Subfolder list starts from the second argument

csv_file="/sci/home/mali.tsadok/pipeline_scripts_singleton.csv"

bash_script="calc_statistic.sh"

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

#subfolder_list=("chr2" "chr3")

# Generate chromosome pairs
pairs=()
for ((i=0; i<${#subfolder_list[@]}; i++)); do
    for ((j=i+1; j<${#subfolder_list[@]}; j++)); do
        pairs+=("${subfolder_list[i]}" "${subfolder_list[j]}")
    done
done

num_pairs=$((${#pairs[@]} / 2))


# Create a temporary SLURM array script
cat <<- EOF > /tmp/sbatch_array.sh
#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --mem=80G

pairs=(${pairs[@]})  # Use the passed subfolder list

chr1=\${pairs[\$((SLURM_ARRAY_TASK_ID-1))*2]}  # Use the SLURM_ARRAY_TASK_ID to get the current subfolder
chr2=\${pairs[\$((SLURM_ARRAY_TASK_ID-1))*2 +1]}  # Use the SLURM_ARRAY_TASK_ID to get the current subfolder

# Run the Python script for the current subfolder
echo "Processing subfolder: \$chr1"
echo "Processing subfolder: \$chr2"

echo "start"
echo $script 
echo "$output_dir_full"
python /sci/labs/orzuk/mali.tsadok/UKB/singleton_code/$script "\$chr1" "\$chr2" "$input_dir_full" "$output_dir_full"

echo "end"
EOF


# Submit the SLURM array job
sbatch --array=1-$num_pairs /tmp/sbatch_array.sh "${pairs[@]}"

echo "Array job submitted with $num_subfolders tasks!"
