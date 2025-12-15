import os 
import pandas as pd
import numpy as np
import sys 
from biomart import BiomartServer
from io import StringIO
import requests
# def ensembl (gene_list, filepath):
    
#     ############# try this pip install pyensembl
    
#     ###find another package or replace to R .

    
    
#     print ("initial server")
#     server = BiomartServer("http://www.ensembl.org/biomart")
#     server.verbose = False
    
#     dataset='hsapiens_gene_ensembl'
#     ensembl = server.datasets[dataset]
    
    
#     response  = ensembl.search({'attributes': ['ensembl_gene_id', 'external_gene_name' ,'chromosome_name' , 'start_position' , 'end_position' ], 
#        'filters': {'ensembl_gene_id':gene_list}})
#     output_string = response.raw.data.decode('ascii')
    
#     output_string_with_headers = "gene\tsymbol\tchromosome\tstart\tend\n" + output_string
      
#     fake_file = StringIO(output_string_with_headers)
#     ensembel_df = pd.read_csv(fake_file, delimiter='\t')
#     ensembel_df['location'] = ensembel_df['chromosome'].astype(str) + ':' + ensembel_df['start'].astype(str) + '-' + ensembel_df['end'].astype(str)
#     return ensembel_df 



def fetch_ensembl_data(gene_list, batch_size=100):
    url = "https://rest.ensembl.org/lookup/id"
    headers = {"Content-Type": "application/json"}

    all_results = []

    # Split the list into batches
    for i in range(0, len(gene_list), batch_size):
        batch = gene_list[i:i + batch_size]

        response = requests.post(url, headers=headers, json={"ids": batch})
        
        if response.status_code != 200:
            print(f"Error: {response.status_code} - {response.text}")
            continue  # Skip this batch

        gene_data = response.json()
        
        for gene_id, details in gene_data.items():
            if details:  # Ensure data exists
                all_results.append({
                    "gene": gene_id,
                    "symbol": details.get("display_name", ""),
                    "chromosome": details.get("seq_region_name", ""),
                    "start": details.get("start", ""),
                    "end": details.get("end", ""),
                })

    # Convert results to a DataFrame
    df_ensembl = pd.DataFrame(all_results)
    df_ensembl["location"] = df_ensembl["chromosome"].astype(str) + ":" + df_ensembl["start"].astype(str) + "-" + df_ensembl["end"].astype(str)
    return df_ensembl 
   

def main (chromosome_name ,folder_input ,folder_output):
    for filename in os.listdir(folder_input):
        print (filename)
        if 'dict' in filename:  
            filepath =folder_input +filename
            df = pd.read_csv(filepath)
            gene_list = df.gene.unique().tolist()
            ensembl_df =  fetch_ensembl_data (gene_list)
            df = df.merge(ensembl_df ,left_on = 'gene' , right_on = 'gene', how = 'left', suffixes=('_1', '_2'))
            df.to_csv(folder_output+filename, index=False)
               

if __name__ == "__main__":

    chromosome_name = sys.argv[1]

    print("Parameter 1:" +str(chromosome_name))
    folder_input =  sys.argv[2]+chromosome_name+"/"
    folder_output = sys.argv[3]+chromosome_name+"/"
    print ("folder_output")
    print (folder_output)
    if 'syn' in folder_input :
        type = 2
    else:
         type = 1
    print ("type")
    print(type)
    main (chromosome_name ,folder_input ,folder_output ) 
    
     
    
    