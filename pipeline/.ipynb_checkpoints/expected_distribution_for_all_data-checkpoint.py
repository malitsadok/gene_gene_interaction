# import pandas as pd
# import matplotlib.pyplot as plt
# import os

# if __name__ == "__main__":
#     # List of files and expected columns
#     files_info = [
#         ("/sci/labs/orzuk/mali.tsadok/UKB/all_data/results/final/lof/lof_final_result.csv", "expected_both_lof"),
#         ("/sci/labs/orzuk/mali.tsadok/UKB/all_data/results/final/lof/missense_final_result.csv", "expected_both_missense"),
#         ("/sci/labs/orzuk/mali.tsadok/UKB/all_data/results/final/lof/missense_lof_final_result.csv", "expected_both_missense_lof"),
#         ("/sci/labs/orzuk/mali.tsadok/UKB/all_data/results/final/lof/lof_missense_final_result.csv", "expected_both_lof_missense")
#     ]

#     output_dir = "/sci/labs/orzuk/mali.tsadok/UKB/all_data/"

#     for file_path, expected_col in files_info:
#         try:
#             print(f"Processing {file_path}...")
#             df = pd.read_csv(file_path, usecols=[expected_col])
#             df = df.dropna()

#             # Filter values up to and including 20
#             df = df[df[expected_col] <= 20]

#             plt.figure(figsize=(8, 5))
#             plt.hist(df[expected_col], bins=50, color='skyblue', edgecolor='black')
#             plt.title(f"Histogram of {expected_col} (≤ 20)")
#             plt.xlabel(expected_col)
#             plt.ylabel("Count")

#             base_name = os.path.basename(file_path).replace(".csv", "")
#             output_img = os.path.join(output_dir, f"expected_distribution_{base_name}_matplotlib_upto20.png")
#             plt.tight_layout()
#             plt.savefig(output_img)
#             plt.close()

#             print(f"Saved: {output_img}")

#         except Exception as e:
#             print(f"Error processing {file_path}: {e}")
# import pandas as pd
# import matplotlib.pyplot as plt
# import os

# if __name__ == "__main__":
#     # List of files and expected columns
#     files_info = [
#         ("/sci/labs/orzuk/mali.tsadok/UKB/all_data/results/final/lof/lof_final_result.csv", "expected_both_lof"),
#         ("/sci/labs/orzuk/mali.tsadok/UKB/all_data/results/final/lof/missense_final_result.csv", "expected_both_missense"),
#         ("/sci/labs/orzuk/mali.tsadok/UKB/all_data/results/final/lof/missense_lof_final_result.csv", "expected_both_missense_lof"),
#         ("/sci/labs/orzuk/mali.tsadok/UKB/all_data/results/final/lof/lof_missense_final_result.csv", "expected_both_lof_missense")
#     ]

#     output_dir = "/sci/labs/orzuk/mali.tsadok/UKB/all_data/"

#     for file_path, expected_col in files_info:
#         try:
#             print(f"Processing {file_path}...")
#             df = pd.read_csv(file_path, usecols=[expected_col])
#             df = df.dropna()

#             # Filter values up to and including 20
#             df = df[df[expected_col] <= 20]

#             # Calculate stats
#             mean_val = df[expected_col].mean()
#             median_val = df[expected_col].median()

#             # Plot
#             plt.figure(figsize=(8, 5))
#             plt.hist(df[expected_col], bins=50, color='skyblue', edgecolor='black')
#             plt.yscale('log')  # log scale for better visibility
#             plt.xlim(0, 20)     # explicit x-axis limit

#             # Title includes explanation
#             plt.title(
#                 f"Histogram of {expected_col} (≤ 20 carriers, log Y)\n"
#                 f"Mean = {mean_val:.2f}, Median = {median_val:.2f}\n"
#                 f"Expected number of individuals carrying both variants (gene pair level)"
#             )

#             plt.xlabel("Expected number of individuals")
#             plt.ylabel("Count (log scale)")

#             # Save
#             base_name = os.path.basename(file_path).replace(".csv", "")
#             output_img = os.path.join(output_dir, f"expected_distribution_{base_name}_logscale_upto20.png")
#             plt.tight_layout()
#             plt.savefig(output_img)
#             plt.close()

#             print(f"Saved: {output_img}")

#         except Exception as e:
#             print(f"Error processing {file_path}: {e}")
# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns
# import numpy as np

# # Set Seaborn style
# sns.set(style="whitegrid")

# # Define input files and plot titles
# files_info = [
#     ("/sci/labs/orzuk/mali.tsadok/UKB/all_data/results/final/lof/lof_final_result.csv", "expected_both_lof", "LoF"),
#     ("/sci/labs/orzuk/mali.tsadok/UKB/all_data/results/final/lof/missense_final_result.csv", "expected_both_missense", "Missense"),
#     ("/sci/labs/orzuk/mali.tsadok/UKB/all_data/results/final/lof/missense_lof_final_result.csv", "expected_both_missense_lof", "Missense + LoF"),
#     ("/sci/labs/orzuk/mali.tsadok/UKB/all_data/results/final/lof/lof_missense_final_result.csv", "expected_both_lof_missense", "LoF + Missense"),
# ]

# # Create figure with subplots
# fig, axes = plt.subplots(2, 2, figsize=(14, 10))
# fig.suptitle("Expected Number of Individuals with Both Variants per Gene Pair\n(MAF ≤ 0.1%, Het and Hom Counted Equally)", fontsize=16)

# for ax, (file_path, column, base_title) in zip(axes.flat, files_info):
#     try:
#         df = pd.read_csv(file_path, usecols=[column])
#         df = df.dropna()
#         df = df[df[column] <= 30]  # Limit x-axis to 30

#         values = df[column]
#         mean_val = values.mean()
#         median_val = values.median()

#         # Histogram + KDE
#         sns.histplot(values, bins=50, kde=True, ax=ax,
#                      color="skyblue", edgecolor="black",
#                      log_scale=(False, True))

#         ax.set_xlim(0, 30)
#         ax.set_xlabel("Expected Individuals")
#         ax.set_ylabel("Count (log scale)")
#         ax.set_title(f"{base_title}\nMean: {mean_val:.2f}, Median: {median_val:.2f}")

#     except Exception as e:
#         ax.set_title(f"Error in {base_title}")
#         print(f"Error with {base_title}: {e}")

# # Add explanation as caption below the plot
# caption = (
#     "Each subplot shows the distribution of the expected number of individuals carrying both variants in a gene pair.\n"
#     "The Minor Allele Frequency (MAF) cutoff is ≤ 0.1%.\n"
#     "Both heterozygous and homozygous mutations are counted as one — i.e., if an individual is homozygous, they are still counted only once."
# )
# fig.text(0.5, -0.02, caption, ha='center', va='center', fontsize=11)

# plt.tight_layout(rect=[0, 0.05, 1, 0.93])  # leave space for title and caption

# # Save the figure
# output_path = "/sci/labs/orzuk/mali.tsadok/UKB/all_data/expected_distribution_final_thesis_plot.png"
# plt.savefig(output_path, bbox_inches='tight', dpi=300)
# plt.close()

# print(f"Saved thesis-ready plot to: {output_path}")


import pandas as pd
import matplotlib.pyplot as plt
import os

if __name__ == "__main__":
    # List of files and expected columns
    files_info = [
        ("/sci/labs/orzuk/mali.tsadok/UKB/all_data/results/final/lof/lof_final_result.csv", "expected_both_lof", "LoF"),
        ("/sci/labs/orzuk/mali.tsadok/UKB/all_data/results/final/lof/missense_final_result.csv", "expected_both_missense", "Missense"),
        ("/sci/labs/orzuk/mali.tsadok/UKB/all_data/results/final/lof/missense_lof_final_result.csv", "expected_both_missense_lof", "Missense→LoF"),
        ("/sci/labs/orzuk/mali.tsadok/UKB/all_data/results/final/lof/lof_missense_final_result.csv", "expected_both_lof_missense", "LoF→Missense")
    ]

    output_dir = "/sci/labs/orzuk/mali.tsadok/UKB/all_data/"
    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    axes = axes.flatten()  # Flatten to iterate easily

    for idx, (file_path, expected_col, title) in enumerate(files_info):
        try:
            print(f"Processing {file_path}...")
            df = pd.read_csv(file_path, usecols=[expected_col])
            df = df.dropna()
            df = df[df[expected_col] <= 100]

            mean_val = df[expected_col].mean()
            median_val = df[expected_col].median()

            ax = axes[idx]
            ax.hist(df[expected_col], bins=50, color='skyblue', edgecolor='black')
            ax.set_yscale('log')
            ax.set_xlim(0, 100)
            ax.set_title(
                f"{title}\nMean = {mean_val:.2f}, Median = {median_val:.2f}",
                fontsize=10
            )
            ax.set_xlabel("Expected # of individuals")
            ax.set_ylabel("Count (log scale)")

        except Exception as e:
            print(f"Error processing {file_path}: {e}")

    # Caption under all subplots
    caption = (
        "Each subplot shows the distribution of the expected number of individuals carrying both variants in a gene pair.\n"
        "The Minor Allele Frequency (MAF) cutoff is ≤ 0.1%.\n"
        "Both heterozygous and homozygous mutations are counted as one — i.e., if an individual is homozygous, they are still counted only once."
    )

    # Adjust layout and add caption
    plt.tight_layout(rect=[0, 0.1, 1, 0.95])
    fig.suptitle("Expected Co-carriers of Gene Pairs (≤ 100)", fontsize=14)
    fig.text(0.5, 0.01, caption, ha='center', va='bottom', fontsize=10)

    # Save combined figure
    output_img = os.path.join(output_dir, "combined_expected_distribution_logscale_upto100.png")
    plt.savefig(output_img)
    plt.close()
    print(f"Saved: {output_img}")

