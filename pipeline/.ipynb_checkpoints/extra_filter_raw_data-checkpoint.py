import os
import sys
import pandas as pd

if __name__ == "__main__":

    chromosome_name = sys.argv[1]
    base_path = "/sci/labs/orzuk/mali.tsadok/UKB/raw_data/lof"
    print ("chromosome_name") 
    print (chromosome_name)
    path = base_path +"/" + chromosome_name
    print ("path")
    print (path)
    all_files = os.listdir(path)
    dict_files = [file for file in all_files if "dict" in file]
    data_files = [file for file in all_files if "data" in file]

    # Mutation types of interest
    options = [
        "missense_variant&splice_region_variant",
        "missense_variant&disruptive_inframe_insertion",
        "missense_variant&conservative_inframe_insertion"
    ]

    # Load only the needed columns from Alpha to save memory
    alpha = pd.read_csv(
        "/sci/labs/orzuk/mali.tsadok/UKB/AlphaMissense_likely_pathogenic.csv",
        usecols=["SNV", "transcript_id"]
    )

    for dict_file in dict_files:
        # Find the matching data file
      
        dict_path = os.path.join(path, dict_file)
        dict_df = pd.read_csv(dict_path)

        # Split dict into missense and non-missense
        dict_missense = dict_df[dict_df["mutation_type"].isin(options)]
        dict_other = dict_df[~dict_df["mutation_type"].isin(options)]

        # Save counts before filtering
        missense_before = len(dict_missense)

        # Inner merge missense with alpha
        dict_missense_filtered = dict_missense.merge(
            alpha,
            left_on=["SNV_separated", "transcript"],
            right_on=["SNV", "transcript_id"],
            how="inner"
        )

        missense_after = len(dict_missense_filtered)

        # Concat filtered missense back with all non-missense
        dict_final = pd.concat([dict_other, dict_missense_filtered], ignore_index=True)

        # Save updated dict
        dict_final.to_csv(dict_path, index=False)

       
        # Print summary
        print(f"File: {dict_file}")
        print(f"  Missense rows before: {missense_before}")
        print(f"  Missense rows after : {missense_after}")
        print(f"  Rows dropped        : {missense_before - missense_after}")
        print(f"  Final dict size     : {len(dict_final)}")

        print("-" * 50)
