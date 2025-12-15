import  os 
import pandas as pd
import re
from io import StringIO
import numpy as  np
from scipy.stats import poisson
import itertools 
import sys


def count_appearance_for_each_snv(df) : 

    #Count the number of individuals who have  SNV in their DNA.
    snv_count_df = pd.DataFrame({"count" : df.groupby(['SNV','count_Mut'  ,'gene'  ]).count_Mut.count()} )
    snv_count_df = snv_count_df.pivot_table( index=['SNV', 'gene'], columns='count_Mut', values='count', fill_value=0).reset_index()
    
    snv_count_df.rename(columns={1: 'countHetero' ,2 : "countHomo"}, inplace=True)
    return snv_count_df

def concat_chromosome_files_with_parts(dfs, chromosome_name, folder_output):
    create_folder_if_not_exists(folder_output)
    
    # Sort file keys based on block number (assuming format like 'ukb23157_c21_b2_v1_dict.csv')
    sorted_file_keys = sorted(dfs.keys(), key=lambda x: int(x.split('_b')[1].split('_')[0]))
    
    combined_df = pd.DataFrame()
    all_genes_in_combined = set()
    common_genes = set()

    j = 0 
  
    for i in range(0, len(sorted_file_keys), 10):
   

        current_files = sorted_file_keys[i:i+10] 

        combined_block_df = pd.DataFrame()
        
        # Step 1: Combine data from current blocks
        for file in current_files:
            current_df = dfs[file]
            combined_block_df = pd.concat([combined_block_df, current_df], ignore_index=True)
        if  common_genes : #remove the ones I added to the prevous block  
            combined_block_df = combined_block_df[~combined_block_df['gene'].isin(common_genes)]
        all_genes_in_combined = set(combined_block_df['gene'].unique())
        
        if i + 10 < len(sorted_file_keys):  # If there is a next block
            next_block_file = sorted_file_keys[i + 10]
            next_block_df = dfs[next_block_file]
            
            next_block_genes = set(next_block_df['gene'].unique())

            
            common_genes = all_genes_in_combined.intersection(next_block_genes)
  
            
            if common_genes:
                gene_overlap_rows = next_block_df[next_block_df['gene'].isin(common_genes)]
                print ( gene_overlap_rows)   
                combined_block_df = pd.concat([combined_block_df, gene_overlap_rows], ignore_index=True)
        
        # Step 3: Save the combined data to a CSV
        block_start = i
        block_end = min(i + 9, len(sorted_file_keys) - 1)  # Adjust the block end if fewer than 10 blocks remain
      
        
        print (folder_output + chromosome_name+"_"+str(j)+".csv")
        combined_block_df.to_csv(folder_output+chromosome_name+"_"+str(j)+".csv", index=False)
        j= j+1 
      
        
        # Step 4: Update combined genes
        all_genes_in_combined.update(combined_block_df['gene'].unique())
        
        # Clear the combined DataFrame for the next set of blocks
        combined_df = pd.DataFrame()
    
    print("Concatenation and merging completed.")





def def_count(value) :
    pattern = r'^([0-9.]+)/([0-9.]+)$'
    match = re.match(pattern, value)
    if match:
        if match.group(1) == '.' : 
           num1 = 0 
        else :
            num1 = match.group(1)
        if match.group(2) == '.' : 
            num2 = 0 
        else :
            num2 = match.group(2)
            
        num1  = int (num1)
        num2  = int (num2)
        if num1 > 0 and  num2 > 0 and num1 == num2:
            return 2 
        if num1 > 0 or  num2 > 0 :
            return 1 
    
    return 0



def process_file (name , folder_input , type=1  ) : 
## type 1 for missnece and 2 for syn

    total_individuals_number  = 460692
    
    dtype_dict = {'id': str}
    #participants_uk = pd.read_csv("C:/Users/maliz/thesa/UKbiobank/data/participants_uk.csv",  dtype=dtype_dict)
    participants_uk = pd.read_csv("/sci/labs/orzuk/mali.tsadok/UKB/participants_uk.csv",  dtype=dtype_dict)

 
    all_files = os.listdir(folder_input)
    dict_files =  [file for file in all_files if "dict" in file ]
    data_files = [file for file in all_files if "data" in file]
    dfs = {}
    for dict_file in dict_files :
        
        print (folder_input)
        print (folder_input+dict_file)
        data_file_name = dict_file.replace("dict", "data")

        data_file = [file for file in data_files if data_file_name  in file]
        data_df = pd.read_csv(folder_input+data_file[0] )
        full_dict_df = pd.read_csv(folder_input+dict_file )


        if type == 1 : 
                options = ["frameshift_variant&splice_region_variant" ,"frameshift_variant" , "stop_gained" ,     "frameshift_variant&stop_lost" , "stop_lost" , "stop_gained&splice_region_variant" , "stop_gained&conservative_inframe_insertion" ,"frameshift_variant&stop_gained" , "stop_lost&disruptive_inframe_deletion" , "rare_amino_acid_variant" , "stop_gained&disruptive_inframe_deletion" , "frameshift_variant&stop_gained&splice_region_variant" , "frameshift_variant&start_lost", "stop_lost&conservative_inframe_deletion" , "start_lost" , "frameshift_variant&start_lost&splice_region_variant" ,
                   "frameshift_variant&stop_lost&splice_region_variant" , "start_lost&conservative_inframe_deletion" , "start_lost&conservative_inframe_insertion" ,"start_lost&disruptive_inframe_deletion" "start_lost&disruptive_inframe_insertion" , "start_lost&splice_region_variant" , "stop_gained&disruptive_inframe_deletion&splice_region_variant" ,
"stop_gained&disruptive_inframe_insertion" , "stop_gained&disruptive_inframe_insertion&splice_region_variant" , 
"stop_lost&conservative_inframe_deletion&splice_region_variant" , "stop_lost&disruptive_inframe_insertion" , "stop_lost&splice_region_variant" ,"bidirectional_gene_fusion" , "splice_donor_variant" ,"splice_acceptor_variant" ,"gene_fusion" ,"frameshift_variant&missense_variant" , "frameshift_variant&synonymous_variant" , "stop_gained&conservative_inframe_insertion&splice_region_variant" ,"start_lost&conservative_inframe_deletion&splice_region_variant" ,"stop_lost&disruptive_inframe_deletion&splice_region_variant" ,
                  "stop_gained&conservative_inframe_deletion" ,"start_lost&disruptive_inframe_deletion&splice_region_variant" ,"exon_loss_variant" ,
                  "frameshift_variant&missense_variant&splice_region_variant" , "missense_variant" , "missense_variant&splice_region_variant" , "missense_variant&disruptive_inframe_insertion" ,"missense_variant&conservative_inframe_insertion"  ]

        else : 
            options = ['synonymous_variant' ]

        full_dict_df = full_dict_df[full_dict_df.mutation_type.isin(options)]
        full_dict_df = full_dict_df.drop(columns=["SNV"], errors="ignore")

        full_dict_df.rename(columns={'SNV_separated': 'SNV'} , inplace = True)
        data_df.rename(columns={'SNV_separated': 'SNV'} , inplace = True)
        print (full_dict_df.head())
        count_unique_gene_df  = full_dict_df.groupby(["SNV"]).gene.nunique()
        
        
        print ("number if unique snv is: "+str (count_unique_gene_df.shape[0]) )
        
        more_than_one_gene_df = count_unique_gene_df[count_unique_gene_df >1]
        print ("There are " + str(more_than_one_gene_df.shape[0]) + " snv with more than one gene related")
        
          
        full_dict_df = full_dict_df.groupby('SNV').filter(lambda x: len(x['gene'].unique()) ==  1)
        print ("shape after filtering  all the snv related to more than one gene is : " +str (full_dict_df.shape[0]))

        
        #data_df =  data_df.drop(columns=['transcript'])
        #data_df = data_df.drop_duplicates()
    

        data_df['count_Mut'] = data_df['Mut'].apply(def_count) #change the mut to number

        dict_df =full_dict_df[['SNV' , 'gene','mutation_type' ]].copy()
        dict_df = dict_df.drop_duplicates(subset="SNV")
        print ("data_df.shape")
        print (data_df.shape)
        data_df  = dict_df.merge(data_df, on= ['SNV' , 'transcript'] , how = 'inner')
        print (data_df.shape)

        data_df = data_df.drop_duplicates( subset=['SNV' , 'gene' , 'Individual' , 'transcript']).copy()
        print ("df after drop duplicates" + str(data_df.shape))
        
        data_df_cutoff = data_df.copy()
        data_df_cutoff = data_df_cutoff[data_df_cutoff.Individual.isin(participants_uk.id)]
        print ("df after participants_uk number" + str(data_df_cutoff.shape))
        


        cutoff_frequnecy = 0.001 * total_individuals_number  # 0.1% of population

        # Step 1: Count how many individuals each SNV appears in
        snv_counts = data_df_cutoff.groupby('SNV')['Individual'].nunique()

        # Step 2: Identify SNVs below or equal to cutoff
        valid_snvs = snv_counts[snv_counts <= cutoff_frequnecy].index

       # Step 3: Print stats before and after filtering
        total_snv_before = data_df_cutoff['SNV'].nunique()
        data_df_cutoff_with_frequency = data_df_cutoff[data_df_cutoff['SNV'].isin(valid_snvs)]
        total_snv_after = data_df_cutoff_with_frequency['SNV'].nunique()
        filtered_out = total_snv_before - total_snv_after

        print(f"Total SNVs before filtering: {total_snv_before}")
        print(f"Total SNVs after filtering (≤ 0.1% of individuals): {total_snv_after}")
        print(f"Number of SNVs filtered out: {filtered_out}")
        dfs[dict_file] = data_df_cutoff_with_frequency
       
        
    return dfs
        
        #filtered_df_cutoff.to_csv(folder_output+name+ "_"+ str(number)+ ".csv", index= False)
        
def create_folder_if_not_exists(folder_path):
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
      

if __name__ == "__main__":

    chromosome_name = sys.argv[1]
    print("Parameter 1:" +str( chromosome_name))
    folder_input =  sys.argv[2]+chromosome_name+"/"
    folder_output = sys.argv[3]+chromosome_name+"/"

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
    dfs = process_file (chromosome_name,folder_input, type )
    concat_chromosome_files_with_parts(dfs,chromosome_name ,folder_output)
   





