# import pandas as pd
# import matplotlib.pyplot as plt
# import os

# if __name__ == "__main__":
#     files_info = [
#         (
#             "/sci/labs/orzuk/mali.tsadok/UKB/singleton_1/results/final/lof/lof_final_result.csv",
#             "both_lof",
#             "Both genes are LoF"
#         ),
#         (
#             "/sci/labs/orzuk/mali.tsadok/UKB/singleton_1/results/final/lof/missense_final_result.csv",
#             "both_missense",
#             "Both genes are missense"
#         ),
#         (
#             "/sci/labs/orzuk/mali.tsadok/UKB/singleton_1/results/final/lof/lof_missense_combined_final_result.csv",
#             "both_lof_missense_combined",
#             "One gene is LoF and the other is missense"
#         )
#     ]

#     output_dir = "/sci/labs/orzuk/mali.tsadok/UKB/"

#     # Compact horizontal layout: 1 row × 3 columns
#     fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
#     plt.style.use("seaborn-v0_8-whitegrid")

#     for idx, (file_path, expected_col, title) in enumerate(files_info):
#         try:
#             print(f"Processing {file_path}...")
#             df = pd.read_csv(file_path, usecols=[expected_col]).dropna()

#             # Stats for annotation
#             max_val = df[expected_col].max()
#             mean_val = df[expected_col].mean()
#             median_val = df[expected_col].median()
#             # Count how many times each integer 0–5 appears
#             counts_0_5 = df[expected_col].value_counts().reindex(range(0, 6), fill_value=0)
#             print(f"\nValue counts for {expected_col} (0–5):")
#             print(counts_0_5)

#             ax = axes[idx]
#             bins = range(int(df[expected_col].min()), int(df[expected_col].max()) + 2)
#             ax.hist(df[expected_col], bins=bins, color="#DAA520", edgecolor="black", alpha=0.8)

#             max_x = int(df[expected_col].max())
#             ax.set_xticks(range(0, max_x + 1, 5))
#             ax.set_xlim(-0.5, 40.5)

#             #ax.hist(df[expected_col], bins=50, color="#DAA520", edgecolor="black", alpha=0.8)
#             ax.set_yscale("log")
           

#             # Titles + axis labels (smaller font for compact layout)
#             ax.set_title(
#                 f"{title}\nMean = {mean_val:.2f}, Median = {median_val:.2f}",
#                 fontsize=11,
#                 pad=8,
#                 fontweight='bold'
#             )
#             ax.set_xlabel("Observed # of Individuals", fontsize=9)
#             if idx == 0:
#                 ax.set_ylabel("Count (log scale)", fontsize=9)
#             else:
#                 ax.set_ylabel("")

#             # Smaller ticks for clarity in compact plot
#             ax.tick_params(axis='both', labelsize=8)

#         except Exception as e:
#             print(f"Error processing {file_path}: {e}")

#     # Overall title
#     fig.suptitle("Observed Co-carriers of Gene Pairs", fontsize=14, fontweight='bold')

#     # Save high-resolution figure (publication-quality)
#     output_img = os.path.join(output_dir, "observed_distribution_logscale_upto_max_small.png")
#     plt.savefig(output_img, dpi=600, bbox_inches='tight')
#     plt.close()
#     print(f"Saved: {output_img}")



import pandas as pd
import matplotlib.pyplot as plt
import os

if __name__ == "__main__":
    files_info = [
        (
            "/sci/labs/orzuk/mali.tsadok/UKB/all_data/results/final/lof/lof_final_result.csv",
            "expected_both_lof",
            "Both genes are LoF"
        ),
        (
            "/sci/labs/orzuk/mali.tsadok/UKB/all_data/results/final/lof/missense_final_result.csv",
            "expected_both_missense",
            "Both genes are missense"
        ),
        (
            "/sci/labs/orzuk/mali.tsadok/UKB/all_data/results/final/lof/lof_missense_combined_final_result.csv",
            "expected_both_lof_missense_combined",
            "One gene is LoF and the other is missense"
        )
    ]

    output_dir = "/sci/labs/orzuk/mali.tsadok/UKB/"

    # Compact horizontal layout: 1 row × 3 columns
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
    plt.style.use("seaborn-v0_8-whitegrid")

    for idx, (file_path, expected_col, title) in enumerate(files_info):
        try:
            print(f"Processing {file_path}...")
            df = pd.read_csv(file_path, usecols=[expected_col]).dropna()

            # Stats for annotation
            max_val = df[expected_col].max()
            mean_val = df[expected_col].mean()
            median_val = df[expected_col].median()
            # Count how many times each integer 0–5 appears
            counts_0_5 = df[expected_col].value_counts().reindex(range(0, 6), fill_value=0)
            print(f"\nValue counts for {expected_col} (0–5):")
            print(counts_0_5)
           
            ax = axes[idx]
            # bins = range(int(df[expected_col].min()), int(df[expected_col].max()) + 2)
            # ax.hist(df[expected_col], bins=bins, color="#DAA520", edgecolor="black", alpha=0.8)

            max_x = int(df[expected_col].max())
            

            ax.hist(df[expected_col], bins=50, color="#f794db", edgecolor="black", alpha=0.8)
            ax.set_yscale("log")
           

            # Titles + axis labels (smaller font for compact layout)
            ax.set_title(
                f"{title}\nMean = {mean_val:.2f}, Median = {median_val:.2f}",
                fontsize=11,
                pad=8,
                fontweight='bold'
            )
            ax.set_xlabel("Expected # of Individuals", fontsize=9)
            if idx == 0:
                ax.set_ylabel("Count (log scale)", fontsize=9)
            else:
                ax.set_ylabel("")

            # Smaller ticks for clarity in compact plot
            ax.tick_params(axis='both', labelsize=8)

        except Exception as e:
            print(f"Error processing {file_path}: {e}")

    # Overall title
    fig.suptitle("Expected Co-carriers of Gene Pairs", fontsize=14, fontweight='bold')

    # Save high-resolution figure (publication-quality)
    output_img = os.path.join(output_dir, "expected_distribution_logscale_all_data.png")
    plt.savefig(output_img, dpi=600, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_img}")
