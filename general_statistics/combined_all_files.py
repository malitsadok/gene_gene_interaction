import os
import pandas as pd



if __name__ == "__main__":
    # Define the parent folder where the subfolders are located
    parent_folder = "/sci/labs/orzuk/mali.tsadok/UKB/Individual_gene/"  # Adjust this to your folder path
    output_file = "/sci/labs/orzuk/mali.tsadok/UKB/combined.csv"  # Name of the resulting CSV file
    
    # Create an empty list to hold DataFrames
    dataframes = []
    
    # Walk through all subfolders and files
    for root, dirs, files in os.walk(parent_folder):
        for file in files:
            if file.endswith(".csv"):  # Only process CSV files
                print (file)
                file_path = os.path.join(root, file)
                # Read each CSV file
                df = pd.read_csv(file_path)
                # Append to the list of DataFrames
                dataframes.append(df)
    
    # Concatenate all DataFrames into one
    if dataframes:
        combined_df = pd.concat(dataframes, ignore_index=True)
        # Save the combined DataFrame to a CSV file
        combined_df.to_csv(output_file, index=False)
        print(f"Combined CSV file saved as {output_file}")
    else:
        print("No CSV files found in the folder.")


