import  os 
import pandas as pd
import numpy as  np
import sys

import itertools


def calculate_cooccurrence(row, left_col, right_col):
    left = row[left_col]
    right = row[right_col]
    if left >= 1 and right >= 1:
        return 1
    else:
        return 0



def gene_gene_interactions_different_chromosomes (df1 ,df2, flag_same_file = False ) : 

  
    gene_gene_observed_df = pd.merge(df1, df2, on='Individual' , suffixes=('_1', '_2'))
    

    gene_gene_observed_df['both_lof_singleton_1'] = gene_gene_observed_df.apply(
    lambda row: calculate_cooccurrence(row, 'lof_singleton_1', 'lof_2'), axis=1)

    gene_gene_observed_df['both_lof_singleton_2'] = gene_gene_observed_df.apply(
    lambda row: calculate_cooccurrence(row, 'lof_1', 'lof_singleton_2'), axis=1)

    gene_gene_observed_df['both_missense_singleton_1'] = gene_gene_observed_df.apply(
    lambda row: calculate_cooccurrence(row, 'missense_singleton_1', 'missense_2'), axis=1)

    gene_gene_observed_df['both_missense_singleton_2'] = gene_gene_observed_df.apply(
    lambda row: calculate_cooccurrence(row, 'missense_1', 'missense_singleton_2'), axis=1)

    gene_gene_observed_df['both_lof_missense_singleton_1'] = gene_gene_observed_df.apply(
    lambda row: calculate_cooccurrence(row, 'lof_singleton_1', 'missense_2'), axis=1)

    gene_gene_observed_df['both_lof_missense_singleton_2'] = gene_gene_observed_df.apply(
    lambda row: calculate_cooccurrence(row, 'lof_1', 'missense_singleton_2'), axis=1)

    gene_gene_observed_df['both_missense_lof_singleton_1'] = gene_gene_observed_df.apply(
    lambda row: calculate_cooccurrence(row, 'missense_singleton_1', 'lof_2'), axis=1)

    gene_gene_observed_df['both_missense_lof_singleton_2'] = gene_gene_observed_df.apply(
    lambda row: calculate_cooccurrence(row, 'missense_1', 'lof_singleton_2'), axis=1)

    gene_gene_observed_df=gene_gene_observed_df.fillna(" ")
    gene_gene_grouped = gene_gene_observed_df.groupby(['gene_1',  'gene_2' ]).agg(
    both_lof_singleton_1=('both_lof_singleton_1', 'sum'),
    both_lof_singleton_2=('both_lof_singleton_2', 'sum'),
    both_missense_singleton_1=('both_missense_singleton_1', 'sum'),
    both_missense_singleton_2=('both_missense_singleton_2', 'sum'),
    both_lof_missense_singleton_1=('both_lof_missense_singleton_1', 'sum'),
    both_lof_missense_singleton_2=('both_lof_missense_singleton_2', 'sum'),
    both_missense_lof_singleton_1=('both_missense_lof_singleton_1', 'sum'),
    both_missense_lof_singleton_2=('both_missense_lof_singleton_2', 'sum')
    
    ).reset_index()

    
    gene_list1 = df1["gene"].unique().tolist()
    gene_list2 = df2["gene"].unique().tolist()
    
    gene_df1 = pd.DataFrame(gene_list1, columns=["gene"])
    gene_df2 = pd.DataFrame(gene_list2, columns=["gene"])

   
    chunk_size = 10000  
    chunks1 = [gene_df1[i:i+chunk_size] for i in range(0, len(gene_df1), chunk_size)]
    chunks2 = [gene_df2[i:i+chunk_size] for i in range(0, len(gene_df2), chunk_size)]
    
    # Perform cross join for each pair of chunks
    all_gene_gene = pd.DataFrame()
    
    for i in  range(len(chunks1)):
        chunk1 = chunks1[i]
        for j in range(len(chunks2)):
            chunk2 = chunks2[j]             
            merged_chunk = pd.merge(chunk1, chunk2, how='cross', suffixes=('_1', '_2'))
            if flag_same_file  :
                 merged_chunk = merged_chunk[merged_chunk["gene_1"] != merged_chunk["gene_2"]]
                 merged_chunk["pair"] = merged_chunk.apply(lambda row: tuple(sorted([row["gene_1"], row["gene_2"]])), axis=1)
                 merged_chunk = merged_chunk.drop_duplicates(subset="pair").drop(columns=["pair"])
            all_gene_gene = pd.concat([all_gene_gene , merged_chunk], ignore_index=True )
            

 

    all_gene_gene =pd.merge(all_gene_gene ,   gene_gene_grouped  , left_on = ['gene_1' ,'gene_2'] , right_on =  ['gene_1' ,'gene_2'], how = 'left', suffixes=('_1', '_2'))    
    print ("all_gene_gene.shape[0]")
    print (all_gene_gene.shape[0])

    all_gene_gene.fillna(0, inplace = True)
    return all_gene_gene


if __name__ == "__main__":
    flag_same_file = False
    param1 = sys.argv[1]
    param2 = sys.argv[2]
   

    input_directory = sys.argv[3]
    output_directory = sys.argv[4]
    print("Parameter 1:", param1)
    print("Parameter 2:", param2)
    
    chr_name1 = param1.split('_')[0]
    chr_name2 = param2.split('_')[0]

    if chr_name1 == chr_name2 : 
        flag_same_file = True

    file1 = os.path.join(output_directory, f"all_gene_gene_{param1}_{param2}.csv")
    file2 = os.path.join(output_directory, f"all_gene_gene_{param2}_{param1}.csv")

    if os.path.exists(file1) or os.path.exists(file2):
        print("At least one of the files exists.")
    else:
        df1 = pd.read_csv(input_directory+chr_name1+"/"+param1+".csv" )
        df2  =pd.read_csv(input_directory+chr_name2+"/"+param2+".csv")
    
        all_gene_gene  = gene_gene_interactions_different_chromosomes (df1 ,df2 , flag_same_file )
        all_gene_gene.to_csv(output_directory+ "all_gene_gene_"+param1+"_"+param2+".csv" , index= False)

