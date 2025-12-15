#!/bin/bash
#SBATCH --time=50:00:00
#SBATCH --ntasks=2
#SBATCH --mem=80G

main_directory="/sci/labs/orzuk/mali.tsadok/UKB/Individual_gene/lof"
results_directory="/sci/labs/orzuk/mali.tsadok/UKB/counted_data/lof"

# Get the list of subfolders
subfolders=($(ls -d "$main_directory"/*/))

# Loop through the subfolders in pairs
for ((i = 0; i < ${#subfolders[@]}; i++)); do
  for ((j = i + 1; j < ${#subfolders[@]}; j++)); do
    chr1=$(basename "${subfolders[$i]}")
    chr2=$(basename "${subfolders[$j]}")
	
    echo $chr1
    echo $chr2

    path1="${subfolders[$i]}"
    path2="${subfolders[$j]}"

    # Populate list1 with all the file names from path1 without the .csv extension
    list1=()
    while IFS= read -r -d '' file; do
      filename=$(basename "$file" .csv)  # Extract filename without the extension
      list1+=("$filename")
    done < <(find "$path1" -type f -name "*.csv" -print0)

    # Populate list2 with all the file names from path2 without the .csv extension
    list2=()
    while IFS= read -r -d '' file; do
      filename=$(basename "$file" .csv)  # Extract filename without the extension
      list2+=("$filename")
    done < <(find "$path2" -type f -name "*.csv" -print0)

    # Check for each pair of files
    for file1 in "${list1[@]}"; do
      for file2 in "${list2[@]}"; do
      
		
		result_file3="${results_directory}/all_gene_gene_${file1}_${file2}.csv"
		result_file4="${results_directory}/all_gene_gene_${file2}_${file1}.csv"
		

		if [ ! -f "$result_file3" ] && [ ! -f "$result_file4" ]; then
          echo "Result file $result_file3 do not   exists. check why." 
        fi
      done
    done

  done
done
