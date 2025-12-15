import  os 
import pandas as pd
import numpy as  np
import sys

import itertools

total_individuals_number  = 460692 #number of pepole in general 


def gene_count (input_directory ,chr_name  ) : 
    directory_path  =  input_directory+chr_name
    files = [f for f in os.listdir(directory_path) if os.path.isfile(os.path.join(directory_path, f))]
    dfs = []
    dfs_individual = [] 
    for file_path in files:
        print (file_path)
        df = pd.read_csv(directory_path+"/"+file_path)
        df_grouped = df.groupby("gene")[["lof" ,"missense" , "lof_or_missense"]].sum().reset_index()
        df_grouped_individual = df.groupby("Individual")[["lof" ,"missense" , "lof_or_missense"]].sum().reset_index()
        dfs.append(df_grouped)
        dfs_individual.append(df_grouped_individual)

    final_gene =  pd.concat(dfs, ignore_index=True)
    final_individual =  pd.concat(dfs_individual, ignore_index=True)
    final_individual = final_individual.groupby("Individual").sum().reset_index()
                       
    return final_gene , final_individual


if __name__ == "__main__":
    chr_name = sys.argv[1]
    input_directory = sys.argv[2]
    output_directory = sys.argv[3]
    print("Parameter 1:", chr_name)
   
    
    gene_count_df, count_individual_df  =  gene_count (input_directory ,chr_name  )
    gene_count_df.to_csv(output_directory+ "gene_count_"+chr_name+".csv" , index= False)
    count_individual_df.to_csv(output_directory+ "Individual_count_"+chr_name+".csv" , index= False)
 
