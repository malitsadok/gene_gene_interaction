import os
import sys
import pandas as pd
import statsmodels.stats.multitest as smm
from scipy import stats


# def process_p_values(df, p_value_columns):
#     """Adjusts p-values using FDR correction and adds q-values."""
#     for col in p_value_columns:
#         if col in df:
#             _, q_values, _, _ = smm.multipletests(df[col], alpha=0.05, method="fdr_bh")
#             df[col.replace("p-value", "q-value")] = q_values
#     return df


# def compute_z_scores(df, p_value_columns):
#     """Computes z-scores from p-values and adds them to the DataFrame."""
#     for col in p_value_columns:
#         if col in df:
#             z_col_name = col.replace("p-value", "z-score")
#             df[z_col_name] = stats.norm.ppf(1 - df[col])
#     return df


def main(input_directory, output_directory, key ,columns_to_load ):
    # Ensure input directory exists
    if not os.path.exists(input_directory):
        raise FileNotFoundError(f"Input directory '{input_directory}' does not exist.")

   

    csv_files = [f for f in os.listdir(input_directory) if f.endswith(".csv")]
    if not csv_files:
        raise ValueError("No CSV files found in the input directory.")

    # Read and merge all CSV files
    df_list = []
    for file in csv_files:
        file_path = os.path.join(input_directory, file)
        tmp = pd.read_csv(file_path, usecols=columns_to_load)
        print (tmp[(tmp.gene_1 == 'ENSG00000154358' ) & (tmp.gene_2 == 'ENSG00000146839' ) ])
        df_list.append(tmp)
 

    if not df_list:
        raise ValueError("No valid CSV files could be read.")

    df_merged = pd.concat(df_list, ignore_index=True)
    df_merged = df_merged.groupby(['gene_1',  'gene_2']).first().reset_index()

    # List of p-value columns
    p_value_columns = [
        "p-value_missense_singleton_1",
        "p-value_missense_singleton_2",
        "p-value_lof_singleton_1",
        "p-value_lof_singleton_2",
        "p-value_missense_lof_singleton_1",
        "p-value_missense_lof_singleton_2",
        "p-value_lof_missense_singleton_1",
        "p-value_lof_missense_singleton_2",
        "p-value_missense" ,
        "p-value_lof" ,
        "p-value_missense_lof" ,
        "p-value_lof_missense" ,
        "p-value_lof_missense_2" , 
    ]

   
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)


    # Save final result
    output_file = os.path.join(output_directory, key+"_final_result.csv")
    df_merged.to_csv(output_file, index=False)
    print (df_merged[(df_merged.gene_1 == 'ENSG00000154358' ) & (df_merged.gene_2 == 'ENSG00000146839' ) ])
    print(f"Results saved to {output_file}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py <input_directory> <output_directory>")
        sys.exit(1)

    input_directory = sys.argv[1]
    output_directory = sys.argv[2]
    dict = {"missense_singleton_1" : ["gene_1" , "gene_2","missense_singleton_1","missense_2" ,"both_missense_singleton_1" ,"expected_both_missense_singleton_1" , "p-value_missense_singleton_1"], 
            "missense_singleton_2" : ["gene_1" , "gene_2" ,"missense_1" ,"missense_singleton_2","both_missense_singleton_2" ,"expected_both_missense_singleton_2" , "p-value_missense_singleton_2"],
            "lof_singleton_1" : ["gene_1" , "gene_2","lof_singleton_1" ,"lof_2" ,"both_lof_singleton_1" ,"expected_both_lof_singleton_1" , "p-value_lof_singleton_1"],
            "lof_singleton_2" : ["gene_1" , "gene_2" ,"lof_1" , "lof_singleton_2" ,"both_lof_singleton_2" ,"expected_both_lof_singleton_2" , "p-value_lof_singleton_2"],
            "missense_lof_singleton_1" : ["gene_1" , "gene_2" ,"missense_singleton_1" ,"lof_2" ,"both_missense_lof_singleton_1" ,"expected_both_missense_lof_singleton_1" , "p-value_missense_lof_singleton_1"], 
            "missense_lof_singleton_2" : ["gene_1" , "gene_2", "missense_1" ,  "lof_singleton_2","both_missense_lof_singleton_2" ,"expected_both_missense_lof_singleton_2" , "p-value_missense_lof_singleton_2"]  ,
            "lof_missense_singleton_1" : ["gene_1" , "gene_2" ,"lof_singleton_1","missense_2" ,"both_lof_missense_singleton_1" ,"expected_both_lof_missense_singleton_1" , "p-value_lof_missense_singleton_1"], 
            "lof_missense_singleton_2" : ["gene_1" , "gene_2" , "lof_1" ,"missense_singleton_2" ,"both_lof_missense_singleton_2" ,"expected_both_lof_missense_singleton_2" , "p-value_lof_missense_singleton_2"] ,
           "missense" : ["gene_1" , "gene_2" , "missense_1" , "missense_singleton_2"   , "missense_singleton_1" ,   "missense_2" , "both_missense" , "expected_both_missense" , "p-value_missense" ]  ,
           "lof" : ["gene_1" , "gene_2" , "lof_singleton_1"  ,"lof_2" , "lof_1",  "lof_singleton_2" ,  "both_lof" , "expected_both_lof" , "p-value_lof" ] ,
           "lof_missense_combined" : ["gene_1" , "gene_2" , "lof_singleton_1"  , "missense_2" , "lof_1" ,   "missense_singleton_2" ,  "missense_singleton_1" ,          "lof_2"  , "missense_1" , "lof_singleton_2" , 'both_lof_missense_combined' , "expected_both_lof_missense_combined" , "p-value_lof_missense_combined" ] ,
           "lof_missense" : ["gene_1" , "gene_2" , "lof_singleton_1"  ,"missense_2" , "lof_1",  "missense_singleton_2" ,  "both_lof_missense" ,           "expected_both_lof_missense" , "p-value_lof_missense" ] ,
           "missense_lof" : ["gene_1" , "gene_2" , "missense_singleton_1"  ,"lof_2" , "missense_1",  "lof_singleton_2" ,  "both_missense_lof" , "expected_both_missense_lof" , "p-value_missense_lof" ] } 
           
    for key in dict :
        main(input_directory, output_directory , key , dict[key])
