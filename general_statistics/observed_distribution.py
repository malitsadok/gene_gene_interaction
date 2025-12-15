import pandas as pd
import matplotlib.pyplot as plt
import os

if __name__ == "__main__":
    # List of files and expected columns
    files_info = [
        ("/sci/labs/orzuk/mali.tsadok/UKB/results/final/lof/lof_singleton_1_final_result.csv", "both_lof_singleton_1", "Gene 1 & Gene 2 LoF"),
        ("/sci/labs/orzuk/mali.tsadok/UKB/results/final/lof/missense_singleton_1_final_result.csv", "both_missense_singleton_1", "Gene 1 & Gene 2 Missense"),
        ("/sci/labs/orzuk/mali.tsadok/UKB/results/final/lof/missense_lof_singleton_1_final_result.csv", "both_missense_lof_singleton_1", "Gene 1 Missense & Gene 2 LoF"),
        ("/sci/labs/orzuk/mali.tsadok/UKB/results/final/lof/lof_missense_singleton_1_final_result.csv", "both_lof_missense_singleton_1", "Gene 1 LoF & Gene 2 Missense")
    ]

    output_dir = "/sci/labs/orzuk/mali.tsadok/UKB/"
    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    axes = axes.flatten()  # Flatten to iterate easily

    # Soft color palette

    for idx, (file_path, col, title) in enumerate(files_info):
        try:
            print(f"Processing {file_path}...")
            df = pd.read_csv(file_path, usecols=[col])
            df = df.dropna()
            print(df.shape)
            df = df[df[col] < 100]
            print(df.shape)

            mean_val = df[col].mean()
            median_val = df[col].median()

            ax = axes[idx]
            ax.hist(df[col], bins=50, color= "#fdbf6f", edgecolor='black')
            ax.set_yscale('log')
            ax.set_xlim(0, 100)
            ax.set_title(
                f"{title}\nMean = {mean_val:.2f}, Median = {median_val:.2f}",
                fontsize=11
            )
            ax.set_xlabel("# of individuals")
            ax.set_ylabel("Count (log scale)")

        except Exception as e:
            print(f"Error processing {file_path}: {e}")

    fig.suptitle("Observed Co-carriers of Gene Pairs", fontsize=16, weight="bold")
    plt.tight_layout()
    fig.subplots_adjust(top=0.88)  # Leave room for suptitle

    # Save combined figure
    output_img = os.path.join(output_dir, "combined_observed_distribution_logscale_upto100.png")
    plt.savefig(output_img, dpi=300)
    plt.close()
    print(f"Saved: {output_img}")
