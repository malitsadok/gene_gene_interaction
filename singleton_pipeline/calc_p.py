import sys 
import pandas as pd
import os

def calc_p(folder_input ,folder_output ):
    dfs =[]
    for file_name in os.listdir(folder_input):
        if file_name.startswith('Individual_count') and file_name.endswith('.csv'):
            file_path = os.path.join(folder_input, file_name)
            df = pd.read_csv(file_path)
            dfs.append(df)
    m = 0
    for file_name in os.listdir(folder_input):
        if file_name.startswith('gene_count') and file_name.endswith('.csv'):
            file_path = os.path.join(folder_input, file_name)
            df = pd.read_csv(file_path)
            print (file_name)
            print (df.shape[0])
            m = m+ df.shape[0]


    combined_df = pd.concat(dfs)

    print ("m")
    print (m)
    result_df = combined_df.groupby('Individual', as_index=False)[["lof" , "missense" , "missense_singleton" , "lof_singleton"]].sum()


    result_df["p_lof"] = result_df["lof"] / m
    result_df["p_missense"] = result_df["missense"] / m
    result_df["p_missense_singleton"] = result_df["missense_singleton"] / m
    result_df["p_lof_singleton"] = result_df["lof_singleton"] / m

    print (result_df.head())
    
    result_df.to_csv(folder_output+"total_Individual_count.csv", index= False)


if __name__ == "__main__":

    folder_input =  sys.argv[1]
    folder_output = sys.argv[2]

    print (folder_input)
    print (folder_output)
    #folder_path = "/sci/labs/orzuk/mali.tsadok/UKB/raw_data/"+chromosome_name

  #  concat_chromosome_files_with_parts (folder_path , chromosome_name   ,"/sci/labs/orzuk/mali.tsadok/UKB/combined_data/"+chromosome_name+"/" )
    if 'syn' in folder_input :
        type = 2
    else:
        type = 1
    print ("type")
    print(type)
    calc_p(folder_input ,folder_output )
   