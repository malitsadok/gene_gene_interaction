import pandas as pd
import numpy as np

def global_shuffle(df_permuted):
    """Perform a full shuffle of gene_2 and chr_2."""
    shuffled = df_permuted[['gene_2', 'chr_2']].sample(frac=1, random_state=None).reset_index(drop=True)
    df_permuted[['gene_2', 'chr_2']] = shuffled
    return df_permuted


def permute_gene2(df, number):
    df_permuted = df.copy()
    
    df_permuted = global_shuffle(df_permuted)
    attempt = 0
    while True:
        
        same_chr_mask = df_permuted['chr_1'] == df_permuted['chr_2']
        n_conflicts = same_chr_mask.sum()
        print(f"Attempt {attempt}: {n_conflicts} conflicts")
    
        if n_conflicts == 0:
            out_path = (
                    f"/sci/labs/orzuk/mali.tsadok/UKB/all_data/permutations/"
                    f"permutation_olida_full_version/premutation_number_{number}.csv"
                )
            df_permuted[['gene_1', 'chr_1', 'gene_2', 'chr_2']].to_csv(out_path, index=False)
            print(f"✅ Saved permutation {number}")
            break
        conflict_idx = df_permuted.index[same_chr_mask]
        conflict_chrs = df_permuted.loc[same_chr_mask, 'chr_1'].unique()
    
        if n_conflicts == 1 or len(conflict_chrs) == 1 or attempt > 1000 :
            df_permuted = global_shuffle(df_permuted)
            continue

        conflict_block = df_permuted.loc[conflict_idx, ['gene_2', 'chr_2']].sample(frac=1).reset_index(drop=True)
        df_permuted.loc[conflict_idx, ['gene_2', 'chr_2']] = conflict_block.values

            
               
if __name__ == "__main__":
    df_paralogs = pd.read_csv(
        "/sci/labs/orzuk/mali.tsadok/UKB/all_data/permutations/olida_pairs_full_version.csv"
    )

    num_permutations = 1000
    for number in range(num_permutations):
        permute_gene2(df_paralogs, number)
        number = number +1 

