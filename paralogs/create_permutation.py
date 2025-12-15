import pandas as pd
import numpy as np

def permute_gene2(df , number):
    df_permuted = df.copy()
    attempts = 0
    while True:
        # Shuffle gene_2 and chr_2 together to keep them matched
        shuffled = df[['gene_2', 'chr_2']].sample(frac=1, random_state=None).reset_index(drop=True)
        
        df_permuted['gene_2'] = shuffled['gene_2']
        df_permuted['chr_2'] = shuffled['chr_2']
        
        # Check for any rows where chr_1 == chr_2
        same_chr = df_permuted['chr_1'] == df_permuted['chr_2']
        if not same_chr.any(): 
            print (attempts)
            df_per = df_permuted[['gene_1' , 'chr_1' , 'gene_2' , 'chr_2']].copy()
            df_per.to_csv("/sci/labs/orzuk/mali.tsadok/UKB/all_data/permutations/permutation_paralogs/premutation_number_"+str(number)+".csv" , index = False)
            break 
        
        attempts += 1
    
    return None


if __name__ == "__main__":

    df_paralogs = pd.read_csv("/sci/labs/orzuk/mali.tsadok/UKB/all_data/results/final/lof/paralogs.csv")

    num_permutations = 10000
    permuted_dfs = []
    number = 0
    for i in range(num_permutations):   
        permute_gene2(df_paralogs ,number )
        number = number +1 