

import  scipy.stats as stats
import  os 
import pandas as pd
import numpy as  npf
import numpy as np
import sys
from scipy.stats import chi2_contingency
from scipy.stats import poisson
from datetime import datetime
import itertools
import re


total_individuals_number  = 431433 #number of pepole in general 
# m = 3053 +1108
m = 17676





def fix_duplicate_couples(name1 , name2, grouped_gene_gene_df , type = 1):
     need_to_check_duplicates = grouped_gene_gene_df[grouped_gene_gene_df.count_dup_pair > 1 ] 
     print ("need_to_check_duplicates")
     print (need_to_check_duplicates)
     return grouped_gene_gene_df


   # Function to check the chromosome order in filenames


 
def concat_result_files(name1, name2 , input_dir, output_dir , type=1): 
    # name = name1+ "_" + name2
    # name3 =  name2+ "_" + name1
    # if os.path.isfile( output_dir +name+".csv") or os.path.isfile( output_dir +name3+".csv")  :
    #     return 

    
    files = os.listdir(input_dir)
    pattern = rf"all_gene_gene_({name1}_\d+_{name2}_\d+|{name2}_\d+_{name1}_\d+)"

    data_gene_gene_files = [f for f in files if re.match(pattern, f)]
   
    data_gene_file = [f for f in files if f.startswith('gene_count')]
    Individual_df =pd.read_csv(input_dir+"total_Individual_count.csv")
    
    dfs = []
    print ("data_gene_gene_files")
    print (data_gene_gene_files)
    
    
    for filename in data_gene_gene_files:
        if name1.split('_')[0] in filename or name2.split('_')[0] in filename:
        #if name1 in filename and  name2 in filename :
            file_path = os.path.join(input_dir, filename)
            df = pd.read_csv(file_path)
            dfs.append(df)
    combined_gene_gene_df = pd.concat(dfs, ignore_index=True)

    dfs = []
    
    for filename in data_gene_file:
        if name1.split('_')[0] in filename or name2.split('_')[0] in filename:
        #if name1 in filename or  name2 in filename :
            file_path = os.path.join(input_dir, filename)
            df = pd.read_csv(file_path)
            dfs.append(df)
    gene_df = pd.concat(dfs, ignore_index=True)

    gene_df = gene_df[['gene' ,'lof' ,'missense'  ]]

     
    combined_gene_gene_df.fillna(' ', inplace = True)

    print ("combined_gene_gene_df shape")
    print(combined_gene_gene_df.shape)




    pairs_counts =pd.merge(combined_gene_gene_df , gene_df  , left_on = ['gene_1'] , right_on =  ['gene'], how = 'left', suffixes=('_1', '_2'))
    pairs_counts.drop(columns = ["gene"] ,inplace=True )
    pairs_counts =pd.merge(pairs_counts , gene_df  , left_on = ['gene_2'] , right_on =  ['gene'], how = 'left', suffixes=('_1', '_2'))
    pairs_counts.drop(columns = ["gene"] ,inplace=True )

    pairs_counts["q_missense_1"] = pairs_counts["missense_1"]/total_individuals_number
    pairs_counts["q_lof_1"] = pairs_counts["lof_1"]/total_individuals_number
    pairs_counts["q_missense_1"] = pairs_counts["missense_1"]/total_individuals_number
    pairs_counts["q_lof_1"] = pairs_counts["lof_1"]/total_individuals_number

    pairs_counts["q_missense_2"] = pairs_counts["missense_2"]/total_individuals_number
    pairs_counts["q_lof_2"] = pairs_counts["lof_2"]/total_individuals_number
    pairs_counts["q_missense_2"] = pairs_counts["missense_2"]/total_individuals_number
    pairs_counts["q_lof_2"] = pairs_counts["lof_2"]/total_individuals_number


    print (datetime.now().time())
    p_missense = np.array(Individual_df["p_missense"])
    p_lof = np.array(Individual_df["p_lof"])


    print (datetime.now().time())
    E_missense = Individual_df["missense"].sum()
    E_lof = Individual_df["lof"].sum()
           

    # Precompute the sums outside (they are constants, not row-dependent)
    sum_missense = np.sum(m * p_missense * (m * p_missense - 1))
    sum_lof = np.sum(m * p_lof * (m * p_lof - 1))
    sum_missense_lof = np.sum(m * p_lof * (m * p_missense - 1))
    
    # Vectorized calculations
    pairs_counts["expected_both_missense"] = (
        pairs_counts["q_missense_1"] * pairs_counts["q_missense_2"] 
        * (total_individuals_number**2 / E_missense**2) 
        * sum_missense
    )
    
    pairs_counts["expected_both_lof"] = (
        pairs_counts["q_lof_1"] * pairs_counts["q_lof_2"] 
        * (total_individuals_number**2 / E_lof**2) 
        * sum_lof
    )
    
    pairs_counts["expected_both_missense_lof"] = (
        pairs_counts["q_missense_1"] * pairs_counts["q_lof_2"] 
        * (total_individuals_number**2 / (E_lof * E_missense)) 
        * sum_missense_lof
    )
    
    pairs_counts["expected_both_lof_missense"] = (
        pairs_counts["q_lof_1"] * pairs_counts["q_missense_2"] 
        * (total_individuals_number**2 / (E_lof * E_missense)) 
        * sum_missense_lof
    )

    pairs_counts["expected_both_lof_missense_combined"] = (
       pairs_counts["expected_both_lof_missense"] +
       pairs_counts["expected_both_missense_lof"]
     
    )
    pairs_counts["both_lof_missense_combined"] = (
       pairs_counts["both_lof_missense"] +
       pairs_counts["both_missense_lof"]
     
    )


    # expected_df = pd.DataFrame(pairs_counts["expected_values"].tolist())
    # pairs_counts = pd.concat([pairs_counts, expected_df], axis=1).drop(columns=["expected_values"])
    pairs_counts['p-value_missense'] = poisson.cdf(pairs_counts['both_missense'], pairs_counts['expected_both_missense'])
    pairs_counts['p-value_lof'] = poisson.cdf(pairs_counts['both_lof'], pairs_counts['expected_both_lof'])
    pairs_counts['p-value_missense_lof'] = poisson.cdf(pairs_counts['both_missense_lof'], pairs_counts['expected_both_missense_lof'])   
    pairs_counts['p-value_lof_missense'] = poisson.cdf(pairs_counts['both_lof_missense'], pairs_counts['expected_both_lof_missense'])
    pairs_counts['p-value_lof_missense_combined'] = poisson.cdf(pairs_counts['both_lof_missense_combined'], pairs_counts['expected_both_lof_missense_combined'])
    
    pairs_counts.to_csv(output_dir +name1+"_"+name2+".csv", index=False)
    return pairs_counts


if __name__ == "__main__":
    chr1 = sys.argv[1]
    chr2 = sys.argv[2]
    input_dir = sys.argv[3]
    output_dir = sys.argv[4]
    type = 1 
    if 'syn' in input_dir : 
        type =2
        

    concat_result_files(chr1,  chr2 ,input_dir,output_dir, type)  