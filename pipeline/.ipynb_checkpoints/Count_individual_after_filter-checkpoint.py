import os
import pandas as pd

# Path to the main folder
main_path = r"/sci/labs/orzuk/mali.tsadok/UKB/filter_data/lof"

# List to store DataFrames
dfs = []

# Walk through all subdirectories
for root, dirs, files in os.walk(main_path):
    for file in files:
        if file.endswith(".csv"):  # change if your files have another extension
            file_path = os.path.join(root, file)
            df = pd.read_csv(file_path)  # load the CSV
            dfs.append(df)

# Concatenate all DataFrames into one big DataFrame
big_df = pd.concat(dfs, ignore_index=True)

# Count unique individuals
unique_count = big_df['Individual'].nunique()

print(f"Total unique individuals: {unique_count}")
