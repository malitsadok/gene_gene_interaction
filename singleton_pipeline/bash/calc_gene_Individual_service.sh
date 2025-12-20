#!/bin/bash
#SBATCH --time=50:00:00
#SBATCH --ntasks=1              
#SBATCH --mem=80G

if [ "$#" -lt 1 ]; then
    echo "Please pass syn or lof"
    exit 1
fi

type=$1
subfolder_list=("${@:2}")  # Subfolder list starts from the second argument

csv_file="/sci/home/mali.tsadok/scripts_singleton_pipeline.csv"

bash_script="calc_gene_Individual_service.sh"

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

# If no subfolders are passed as arguments, get all subfolders in the input directory
if [ "${#subfolder_list[@]}" -eq 0 ]; then
    echo "No subfolder list provided. Using all subfolders in the input directory."
    subfolder_list=($(find "$input_dir_full" -mindepth 1 -maxdepth 1 -type d -exec basename {} \;))
fi

num_subfolders=${#subfolder_list[@]}  # Count how many subfolders were determined


# Create a temporary SLURM array script
cat <<- EOF > /tmp/sbatch_array.sh
#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --mem=80G

subfolder_list=(${subfolder_list[@]})  # Use the passed subfolder list
subdir_name=\${subfolder_list[\$SLURM_ARRAY_TASK_ID-1]}  # Use the SLURM_ARRAY_TASK_ID to get the current subfolder

# Run the Python script for the current subfolder
echo "Processing subfolder: \$subdir_name"
echo "start"
python /sci/labs/orzuk/mali.tsadok/UKB/code/$script "\$subdir_name" "$input_dir_full" "$output_dir_full"
echo "Finished processing subfolder: \$subdir_name"
EOF


# Submit the SLURM array job
sbatch --array=1-$num_subfolders /tmp/sbatch_array.sh "${subfolder_list[@]}"

echo "Array job submitted with $num_subfolders tasks!"
