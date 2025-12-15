import pandas as pd

if __name__ == "__main__":

    # Load CSVs
    df1 = pd.read_csv("/sci/labs/orzuk/mali.tsadok/UKB/singleton_1/results/final/lof/lof_final_result.csv")
    df2 = pd.read_csv("/sci/labs/orzuk/mali.tsadok/UKB/singleton_1/results/final/lof/missense_final_result.csv")
    df3 = pd.read_csv("/sci/labs/orzuk/mali.tsadok/UKB/singleton_1/results/final/lof/lof_missense_combined_final_result.csv")

    # Filter but keep all columns
    df1 = df1[df1["expected_both_lof"] > 5]
    df2 = df2[df2["expected_both_missense"] > 5]
    df3 = df3[df3["expected_both_lof_missense_combined"] > 5]

    # ---- Common in all 3 ----
    common_all = (
        df1.merge(df2, on=["gene_1", "gene_2"], suffixes=("_lof", "_missense"))
           .merge(df3, on=["gene_1", "gene_2"])
    )
    # rename columns from df3 to avoid clashes
    common_all = common_all.add_suffix("_combined").rename(columns={
        "gene_1_combined": "gene_1",
        "gene_2_combined": "gene_2"
    })
    common_all.to_csv(
        "/sci/labs/orzuk/mali.tsadok/UKB/singleton_1/results/final/mutual/mutual.csv",
        index=False
    )

    # ---- Common in each pair ----
    common_12 = df1.merge(df2, on=["gene_1", "gene_2"], suffixes=("_lof", "_missense"))
    common_23 = df2.merge(df3, on=["gene_1", "gene_2"], suffixes=("_missense", "_combined"))
    common_13 = df1.merge(df3, on=["gene_1", "gene_2"], suffixes=("_lof", "_combined"))

    common_12.to_csv(
        "/sci/labs/orzuk/mali.tsadok/UKB/singleton_1/results/final/mutual/mutual_lof_missense.csv",
        index=False
    )
    common_23.to_csv(
        "/sci/labs/orzuk/mali.tsadok/UKB/singleton_1/results/final/mutual/mutual_missense_combined.csv",
        index=False
    )
    common_13.to_csv(
        "/sci/labs/orzuk/mali.tsadok/UKB/singleton_1/results/final/mutual/mutual_lof_combined.csv",
        index=False
    )
