import os
import pandas as pd
from glob import glob
from collections import defaultdict



def process_file(path,  global_uniques, global_singletons , num_genes):
  df = pd.read_csv(path)
  num_genes = df.groupby("gene").nunique() + num_genes
   
  counts = (
        df.groupby("mutation_type")["SNV"]
          .nunique()
          .to_dict()
    )

    # --- Count singleton SNVs per mutation type (appear exactly once in this file) ---
  snv_counts = (
        df.groupby(["mutation_type", "SNV"])
          .size()
          .reset_index(name="count")
    )
  singleton_counts = (
        snv_counts[snv_counts["count"] == 1]
          .groupby("mutation_type")["SNV"]
          .nunique()
          .to_dict()
    )
  for mut_type, cnt in counts.items():
        global_uniques[mut_type] += cnt
  for mut_type, cnt in singleton_counts.items():
        global_singletons[mut_type] += cnt
        if mut_type == 'bidirectional_gene_fusion' : 
           print (snv_counts[snv_counts.mutation_type == mut_type].SNV)

  return global_uniques, global_singletons , num_genes


def process_all_subfolders(root_dir):
    all_files = glob(os.path.join(root_dir, '**/*.csv'), recursive=True)
    global_singletons = defaultdict(int)
    global_uniques = defaultdict(int)
    num_genes = 0 
    for file_path in all_files:
        try:
            global_uniques, global_singletons , num_genes = process_file(
                file_path,  global_uniques, global_singletons , num_genes
            )
        except Exception as e:
            print(f"Failed to process : {e}")

    rows = []  
    print ("num of gene is: ")
    print (num_genes)
    for mut_type in set(global_uniques.keys()) | set(global_singletons.keys()):
        rows.append([
            mut_type,
            global_uniques[mut_type],
            global_singletons[mut_type]
        ])
    
    df_counts = pd.DataFrame(rows, columns=["mutation_type", "unique_snv_count", "singleton_snv_count"])
    df_counts.to_csv("/sci/labs/orzuk/mali.tsadok/UKB/mutation_type_after_filter.csv", index=False)


if __name__ == "__main__":
        process_all_subfolders("/sci/labs/orzuk/mali.tsadok/UKB/filter_data/lof")

# Example usage:
# result = process_all_subfolders("/your/root/folder/path")
# result.to_csv("aggregated_summary.csv", index=False)
