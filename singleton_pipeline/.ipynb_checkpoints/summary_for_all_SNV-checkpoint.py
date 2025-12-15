import os
import pandas as pd
from glob import glob

def process_file(path):
    df = pd.read_csv(path)
    
    # Flag mutation types
    option_missense = ["missense_variant" , "missense_variant&splice_region_variant" , "missense_variant&disruptive_inframe_insertion" ,
                      "missense_variant&conservative_inframe_insertion"  ]
    df['is_missense'] = df['mutation_type'].isin(option_missense)
    df['is_lof'] = ~df['is_missense']
    
    # Group and summarize by SNV
    summary_df = df.groupby('SNV').agg(
        gene=('gene', 'first'),
        mutation_type=('mutation_type' , 'first'),
        n_individuals=('Individual', 'nunique'),
        
        n_missense=('is_missense', 'sum'),
        n_lof=('is_lof', 'sum'),
        
        count_hetrozygous_missense=(
            'count_Mut',
            lambda x: ((x == 1) & df.loc[x.index, 'is_missense']).sum()
        ),
        count_homozygous_missense=(
            'count_Mut',
            lambda x: ((x > 1) & df.loc[x.index, 'is_missense']).sum()
        ),
        
        count_hetrozygous_lof=(
            'count_Mut',
            lambda x: ((x == 1) & df.loc[x.index, 'is_lof']).sum()
        ),
        count_homozygous_lof=(
            'count_Mut',
            lambda x: ((x > 1) & df.loc[x.index, 'is_lof']).sum()
        )
    ).reset_index()
    
    return summary_df

def process_all_subfolders(root_dir):
    all_files = glob(os.path.join(root_dir, '**/*.csv'), recursive=True)
    dfs = []

    for file_path in all_files:
        try:
            df_summary = process_file(file_path)
            df_summary['source_file'] = os.path.basename(file_path)  # optional traceability
            dfs.append(df_summary)
        except Exception as e:
            print(f"Failed to process {file_path}: {e}")
    if dfs:
        tmp =  pd.concat(dfs, ignore_index=True)
        tmp.to_csv("/sci/labs/orzuk/mali.tsadok/UKB/summary_all_SNV")
    else:
        return pd.DataFrame()  # empty fallback


if __name__ == "__main__":
        process_all_subfolders("/sci/labs/orzuk/mali.tsadok/UKB/filter_data/lof")

# Example usage:
# result = process_all_subfolders("/your/root/folder/path")
# result.to_csv("aggregated_summary.csv", index=False)
