import pandas as pd

# Load paralogs
df_paralogs = pd.read_csv("/sci/labs/orzuk/mali.tsadok/UKB/all_data/results/final/lof/paralogs.csv")

# Filtering function
def filter_by_paralogs(df_input, df_paralogs, out_path):
    # Merge on gene_1 + gene_2
    merged1 = df_input.merge(df_paralogs, on=['gene_1', 'gene_2'], how='inner')
    # Merge on reversed pairs
    merged2 = df_input.merge(
      df_paralogs,
      left_on=['gene_1', 'gene_2'],
      right_on=['gene_2', 'gene_1'],
      how='inner'
    )

    # Use the columns from the left side (original input)
    columns_to_keep = [col for col in merged2.columns if col.endswith('_x') or col not in ['gene_1_y', 'gene_2_y']]
    merged2 = merged2[columns_to_keep]

    # Rename _x columns back to original
    merged2 = merged2.rename(columns={'gene_1_x': 'gene_1', 'gene_2_x': 'gene_2'}) 

    filtered_df = pd.concat([merged1, merged2], ignore_index=True)
    filtered_df = filtered_df[df_input.columns]
    # Save
    filtered_df.to_csv(out_path, index=False)

# File paths and output destinations
file_configs = [
    {
        "name": "lof",
        "input_path": "/sci/labs/orzuk/mali.tsadok/UKB/all_data/results/final/lof/lof_final_result.csv",
        "output_path": "/sci/labs/orzuk/mali.tsadok/UKB/all_data/results/final/lof/lof_paralogs_results.csv"
    },
    {
        "name": "missense",
        "input_path": "/sci/labs/orzuk/mali.tsadok/UKB/all_data/results/final/lof/missense_final_result.csv",
        "output_path": "/sci/labs/orzuk/mali.tsadok/UKB/all_data/results/final/lof/missense_paralogs_results.csv"
    },
    {
        "name": "missense_lof",
        "input_path": "/sci/labs/orzuk/mali.tsadok/UKB/all_data/results/final/lof/missense_lof_final_result.csv",
        "output_path": "/sci/labs/orzuk/mali.tsadok/UKB/all_data/results/final/lof/missense_lof_paralogs_results.csv"
    },
    {
        "name": "lof_missense",
        "input_path": "/sci/labs/orzuk/mali.tsadok/UKB/all_data/results/final/lof/lof_missense_final_result.csv",
        "output_path": "/sci/labs/orzuk/mali.tsadok/UKB/all_data/results/final/lof/lof_missense_paralogs_results.csv"
    }
]

# Run the filtering for each file
for config in file_configs:
    print(f"Processing {config['name']}...")
    df = pd.read_csv(config["input_path"])
    filter_by_paralogs(df, df_paralogs, config["output_path"])
    print(f"Saved filtered {config['name']} to {config['output_path']}")
