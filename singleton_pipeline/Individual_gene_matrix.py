import os
import pandas as pd
import sys 



def count_appearance_for_each_snv(df) : 

    #Count the number of individuals who have  SNV in their DNA.
    snv_count_df = pd.DataFrame({"count" : df.groupby(['SNV','count_Mut'  ,'gene','mutation_type'  ]).count_Mut.count()} )
    snv_count_df = snv_count_df.pivot_table( index=['SNV', 'gene','mutation_type'], columns='count_Mut', values='count', fill_value=0).reset_index()
    snv_count_df.rename(columns={1: 'countHetero' ,2 : "countHomo"}, inplace=True)
    return snv_count_df

def process_file(file_path):
    print ("file_path")
    print (file_path)
    df = pd.read_csv(file_path)
    snv_count_df = count_appearance_for_each_snv(df)
    columns = ['Individual', 'gene','count_Mut'  ]
    # Identify singleton variants for missense and lof
    singleton_vars = snv_count_df[(snv_count_df.countHetero == 1) & (snv_count_df.countHomo == 0)]
    singleton_vars_missense = singleton_vars[singleton_vars['mutation_type'].str.contains('missense_variant', na=False)]
    singleton_vars_lof = singleton_vars[~singleton_vars['mutation_type'].str.contains('missense_variant', na=False)]

    # Define conditions for categories

# Define conditions for categories
    categories = {
    'singleton_lof': (df.SNV.isin(singleton_vars_lof.SNV), 0, 1, 0, 0),
    'singleton_missense': (df.SNV.isin(singleton_vars_missense.SNV), 0, 0, 0, 1),
    'missense': (df['mutation_type'].str.contains('missense_variant', na=False), 0, 0, 1, 0),
    'lof': (~df['mutation_type'].str.contains('missense_variant', na=False), 1, 0, 0, 0),
     }

    dfs = []
    for category, (condition, lof, lof_singleton, missense, missense_singleton) in categories.items():
    
        filtered_df = df[condition]
        filtered_no_dup_df = filtered_df.drop_duplicates(subset=['gene', 'Individual'])[columns]
    

        filtered_no_dup_df["lof"] = lof
        filtered_no_dup_df["lof_singleton"] = lof_singleton
        filtered_no_dup_df["missense"] = missense
        filtered_no_dup_df["missense_singleton"] = missense_singleton

        # if category == "missense":
        #     filtered_no_dup_df.loc[filtered_no_dup_df["count_Mut"] == 2, "missense"] = 2

        # if category == "lof":
        #     filtered_no_dup_df.loc[filtered_no_dup_df["count_Mut"] == 2, "lof"] = 2

        dfs.append(filtered_no_dup_df)

 
    return pd.concat(dfs, ignore_index=True)

def process(directory_path , output_directory):
    
    files = [f for f in os.listdir(directory_path) if os.path.isfile(os.path.join(directory_path, f))]
    
    print (files)
    for file_path in files:
        full_path = os.path.join(directory_path, file_path)
        
        result_df = process_file(full_path)
        result_df = result_df.drop("count_Mut", axis=1)  # Remove a column by name
        grouped_df = result_df.groupby(['Individual', 'gene'], as_index=False).max()
        os.makedirs(output_directory, exist_ok=True)
        grouped_df.to_csv(output_directory+"/"+ file_path , index= False)




if __name__ == "__main__":
    chr_name = sys.argv[1]
    input_directory = sys.argv[2]
    output_directory = sys.argv[3]
    print("Parameter 1:", chr_name)
    process( input_directory+chr_name , output_directory+chr_name)   
 
    
 

