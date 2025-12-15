import pandas as pd
import os

def load_file(path, chunk_size=100_000):
    chunks = []
    for chunk in pd.read_csv(path, chunksize=chunk_size):
        chunks.append(chunk)
    return pd.concat(chunks, ignore_index=True)

def load_big_files(base_dir):
    big_files = {
        "lof": "lof_final_result.csv" ,
         "missense": "missense_final_result.csv",
         "missense_lof": "missense_lof_final_result.csv",
         "lof_missense": "missense_lof_final_result.csv",  # adjust if this is a different file!
    }
    ref_dfs = {}

    for key, filename in big_files.items():
        path = os.path.join(base_dir, filename)
        df = load_file(path)
        df["gene_pair"] = df.apply(
            lambda row: tuple(sorted([row["gene_1"], row["gene_2"]])),
            axis=1,
        )
        ref_dfs[key] = df
        print(f"Loaded {key} with shape {df.shape}")

    return ref_dfs



def process_permutation_file(
    perm_path,
    ref_dfs,
):
    # Load permutation file
    df_perm = pd.read_csv(perm_path, usecols=['gene_1','chr_1','gene_2','chr_2'])
    df_perm["gene_pair"] = df_perm.apply(
        lambda row: tuple(sorted([row["gene_1"], row["gene_2"]])),
        axis=1,
    )

    # Start with permutation DataFrame
    df_result = df_perm.copy()

    for key, ref_df in ref_dfs.items():
        # Find relevant pairs
        gene_pairs = set(df_perm["gene_pair"])
        ref_filtered = ref_df[ref_df["gene_pair"].isin(gene_pairs)]

        # Find columns we want to merge
        cols_to_merge = [
            c
            for c in ref_filtered.columns
            if c not in ["gene_1", "gene_2", "gene_pair"]
        ]

        # Drop duplicate columns
        cols_not_already_in_df = [
            c for c in cols_to_merge if c not in df_result.columns
        ]

        if not cols_not_already_in_df:
            # All columns already exist, nothing to merge
            continue

        temp_df = ref_filtered[["gene_pair"] + cols_not_already_in_df]

        # Merge into permutation DF
        df_result = df_result.merge(
            temp_df,
            on="gene_pair",
            how="left",
        )

    print (df_result.head())
    df_result.drop(columns="gene_pair", inplace=True)
    df_result.to_csv(perm_path, index=False)
    



def process_all_permutations(
    perm_dir,
    ref_dfs,
):

    for i in range(5000, 5200):
        target_filename = f"premutation_number_{i}.csv"
        perm_path = os.path.join(perm_dir, target_filename)
        
        if os.path.isfile(perm_path):
            process_permutation_file(
                perm_path,
                ref_dfs,
            )
        else:
            print(f"File not found: {perm_path}")
    # for filename in os.listdir(perm_dir):
    #     if filename.endswith(".csv"):
    #         perm_path = os.path.join(perm_dir, filename)
        

    #         process_permutation_file(
    #             perm_path,
    #             ref_dfs,
    #         )

if __name__ == "__main__":
    base_result_dir = "/sci/labs/orzuk/mali.tsadok/UKB/all_data/results/final/lof"
    permutations_dir = "/sci/labs/orzuk/mali.tsadok/UKB/all_data/results/final/permutation"

    # 1. Load big files ONCE
    ref_dfs = load_big_files(base_result_dir)

    # 2. Process all permutations
    process_all_permutations(
        permutations_dir,
        ref_dfs,
    )
