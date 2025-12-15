import os
import pandas as pd

# Base path where chromosome directories are located
base_path = "/sci/labs/orzuk/mali.tsadok/UKB/raw_data/lof"

# Mutation types to keep
options = [
    "frameshift_variant&splice_region_variant", "frameshift_variant", "stop_gained", "missense_variant",
    "missense_variant&splice_region_variant", "frameshift_variant&stop_lost", "stop_lost",
    "stop_gained&splice_region_variant", "stop_gained&conservative_inframe_insertion",
    "frameshift_variant&stop_gained", "stop_lost&disruptive_inframe_deletion", "rare_amino_acid_variant",
    "stop_gained&disruptive_inframe_deletion", "frameshift_variant&stop_gained&splice_region_variant",
    "frameshift_variant&start_lost", "stop_lost&conservative_inframe_deletion", "start_lost",
    "frameshift_variant&start_lost&splice_region_variant", "frameshift_variant&stop_lost&splice_region_variant",
    "missense_variant&disruptive_inframe_insertion", "start_lost&conservative_inframe_deletion",
    "start_lost&conservative_inframe_insertion", "start_lost&disruptive_inframe_deletion",
    "start_lost&disruptive_inframe_insertion", "start_lost&splice_region_variant",
    "stop_gained&disruptive_inframe_deletion&splice_region_variant", "stop_gained&disruptive_inframe_insertion",
    "stop_gained&disruptive_inframe_insertion&splice_region_variant",
    "stop_lost&conservative_inframe_deletion&splice_region_variant", "stop_lost&disruptive_inframe_insertion",
    "stop_lost&splice_region_variant"
]

# Load participants file (UK subset)
dtype_dict = {'id': str}
participants_uk = pd.read_csv(
    "/sci/labs/orzuk/mali.tsadok/UKB/participants_uk.csv",
    dtype=dtype_dict
)

# Loop over chromosomes 1–22
for chrom in range(21, 23):
    path = os.path.join(base_path, f"chr{chrom}")
    print(f"\n===== Processing chromosome {chrom} =====")

    # -----------------------------
    # Load all "data" CSV files
    # -----------------------------
    data_files = [f for f in os.listdir(path) if "data" in f and f.endswith(".csv")]
    dfs_data = [pd.read_csv(os.path.join(path, f)) for f in data_files]
    big_df = pd.concat(dfs_data, ignore_index=True) if dfs_data else pd.DataFrame()

    # -----------------------------
    # Load all "dict" CSV files
    # -----------------------------
    dict_files = [f for f in os.listdir(path) if "dict" in f and f.endswith(".csv")]
    dfs_dict = [pd.read_csv(os.path.join(path, f)) for f in dict_files]
    dict_df = pd.concat(dfs_dict, ignore_index=True) if dfs_dict else pd.DataFrame()

    if big_df.empty or dict_df.empty:
        print(f"No data or dict files found for chromosome {chrom}")
        continue

    # -----------------------------
    # Merge on SNV_separated column
    # -----------------------------
    merged_df = big_df.merge(dict_df, on="SNV_separated", how="left")

    # -----------------------------
    # Filter rows where mutation_type is in options
    # -----------------------------
    filtered_df = merged_df[merged_df["mutation_type"].isin(options)]

    # -----------------------------
    # Count unique individuals after filtering
    # -----------------------------
    unique_individuals = filtered_df["Individual"].nunique()
    print(f"Unique individuals after filtering: {unique_individuals}")

    # -----------------------------
    # Print mutation types that were filtered out
    # -----------------------------
    filtered_out_types = merged_df.loc[~merged_df["mutation_type"].isin(options), "mutation_type"].unique()
    print("Mutation types that did NOT pass the filter:")
    for mut in sorted(filtered_out_types):
        print(f"  - {mut}")
