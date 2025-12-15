import os
import pandas as pd
import sys 





def process_file(file_path):
    print ("file_path")
    print (file_path)
    df = pd.read_csv(file_path)
    columns = ['Individual', 'gene','count_Mut'  ]
    

    # Define conditions for categories
    option_missense = ["missense_variant" , "missense_variant&splice_region_variant" , "missense_variant&disruptive_inframe_insertion" ,
                      "missense_variant&conservative_inframe_insertion"  ]
    missense_condition = df['mutation_type'].isin(option_missense)
    lof_condition = ~missense_condition
    lof_or_missense_condition = missense_condition | lof_condition  

    categories = {
        'missense': (missense_condition, 0, 1, 0),
        'lof': (lof_condition, 1, 0, 0),
        'lof_or_missense': (lof_or_missense_condition, 0, 0, 1)
    }

    dfs = []
    for category, (condition, lof, missense ,lof_or_missense ) in categories.items():
    
        filtered_df = df[condition].copy()
        filtered_no_dup_df = filtered_df.drop_duplicates(subset=['gene', 'Individual'])[columns]
        filtered_no_dup_df["lof"] = lof
        filtered_no_dup_df["missense"] = missense
        filtered_no_dup_df['lof_or_missense'] = lof_or_missense
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
 
    
 

