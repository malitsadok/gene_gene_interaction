import os
import pandas as pd
from glob import glob


# Function to pick correct row
def pick_row(group):
    # LOF priority
    lof_rows = group[group["mutation_type"].isin(option_lof)]
    if not lof_rows.empty:
        return lof_rows.iloc[0]
    
    # Missense second priority
    missense_rows = group[group["mutation_type"].isin(option_missense)]
    if not missense_rows.empty:
        return missense_rows.iloc[0]
    
    # Otherwise take first
    return group.iloc[0]


def process_joined_file(data_path, dict_path, join_on):
    df_dict = pd.read_csv(dict_path)
    df_data = pd.read_csv(data_path)
    # Join on the given key(s)
    df_data = df_data.drop_duplicates(subset=["SNV_separated", "Individual"], keep="last")
    counts = df_dict["SNV_separated"].value_counts()

    # Separate into groups
    df_two = df_dict[df_dict["SNV_separated"].isin(counts[counts > 1].index)]
    df_one = df_dict[df_dict["SNV_separated"].isin(counts[counts == 1].index)]


    df_two_collapsed = (
       df_two.groupby("SNV_separated", group_keys=False)
          .apply(pick_row)
          .reset_index(drop=True))
    
    df_dict = pd.concat([df_one, df_two_collapsed], ignore_index=True)
    

    # Merge on the given key(s)
    df = pd.merge(df_data, df_dict, on=join_on, how='left')
    



    option_lof = [ 'frameshift_variant', 'frameshift_variant&splice_region_variant', 'frameshift_variant&start_lost',
               'frameshift_variant&stop_gained', 'frameshift_variant&stop_gained&splice_region_variant',
               'frameshift_variant&stop_lost', 'start_lost', 'stop_gained', 'stop_gained&splice_region_variant',
               'stop_lost', 'start_lost&conservative_inframe_deletion',
               'start_lost&conservative_inframe_insertion', 'stop_gained&disruptive_inframe_deletion',
               'stop_lost&conservative_inframe_deletion', 'start_lost&splice_region_variant',
               'stop_gained&disruptive_inframe_deletion&splice_region_variant', 'stop_lost&splice_region_variant',
               'bidirectional_gene_fusion', 'stop_gained&conservative_inframe_insertion',
               'stop_gained&disruptive_inframe_insertion', 'start_lost&disruptive_inframe_deletion',
               'stop_lost&disruptive_inframe_deletion', 'gene_fusion', 'rare_amino_acid_variant',
               'stop_gained&conservative_inframe_insertion&splice_region_variant',
               'stop_gained&disruptive_inframe_insertion&splice_region_variant',
               'frameshift_variant&synonymous_variant',
               'start_lost&conservative_inframe_deletion&splice_region_variant',
               'frameshift_variant&stop_lost&splice_region_variant',
               'start_lost&disruptive_inframe_insertion',
               'frameshift_variant&missense_variant&splice_region_variant',
               'frameshift_variant&start_lost&splice_region_variant',
               'frameshift_variant&missense_variant',
               'stop_lost&conservative_inframe_deletion&splice_region_variant',
               'stop_gained&conservative_inframe_deletion',
               'start_lost&disruptive_inframe_deletion&splice_region_variant',
               'stop_lost&disruptive_inframe_deletion&splice_region_variant',
               'start_lost&conservative_inframe_insertion&splice_region_variant',
               'stop_lost&disruptive_inframe_insertion' ]

    option_missense = ["missense_variant", "missense_variant&splice_region_variant",
                   "missense_variant&disruptive_inframe_insertion",
                   "missense_variant&conservative_inframe_insertion"]

    df = df[(df.mutation_type.isin(options_lof) )|  (df.mutation_type.isin(option_missense)) ]
    
    # Flag mutation types (your existing logic)
    df['is_missense'] = df['mutation_type'].isin(option_missense)
    df['is_lof'] = ~df['is_missense']
    
    # Group and summarize by SNV
    summary_df = df.groupby('SNV_separated').agg(
        gene=('gene', 'first'),
        mutation_type=('mutation_type' , 'first'),
        n_individuals=('Individual', 'nunique'),
        n_missense=('is_missense', 'sum'),
        n_lof=('is_lof', 'sum'),
    ).reset_index()
    
    return summary_df

def process_all_pairs(root_dir, join_on):
    # Find all _data.csv files
    data_files = glob(os.path.join(root_dir, '**/*_data.csv'), recursive=True)
    dfs = []

    for data_path in data_files:
        base_name = os.path.basename(data_path).replace('_data.csv', '')
        dict_path = os.path.join(os.path.dirname(data_path), base_name + '_dict.csv')

        if not os.path.exists(dict_path):
            print(f"Matching dict file not found for {data_path}")
            continue
        
        try:
            df_summary = process_joined_file(data_path, dict_path, join_on)
            df_summary['source_file'] = base_name
            dfs.append(df_summary)
        except Exception as e:
            print(f"Failed to process {base_name}: {e}")

    if dfs:
        combined_df = pd.concat(dfs, ignore_index=True)
        combined_df.to_csv("/sci/labs/orzuk/mali.tsadok/UKB/summary_all_SNV_raw_data.csv", index=False)
        return combined_df
    else:
        return pd.DataFrame()  # empty fallback

if __name__ == "__main__":
    # replace 'SNV' with the actual column(s) you want to join on (can be list)
    process_all_pairs("/sci/labs/orzuk/mali.tsadok/UKB/raw_data/lof", join_on='SNV_separated')
