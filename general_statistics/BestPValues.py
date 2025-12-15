import pandas as pd
import statsmodels.stats.multitest as smm

def load_file(path):
    chunk_size = 100000
    chunks = []
    for chunk in pd.read_csv(path, chunksize=chunk_size):
        chunks.append(chunk)
    df = pd.concat(chunks, ignore_index=True)
    return df 

def process_p_values(df, p_value_columns):
    """Adjusts p-values using FDR correction and adds q-values."""
    for col in p_value_columns:
        if col in df:
            _, q_values, _, _ = smm.multipletests(df[col], alpha=0.05, method="fdr_bh")
            df[col.replace("p-value", "q-value")] = q_values
    return df

if __name__ == "__main__":
    # Define the files to process
    configs = [
        {"file": "missense_final_result.csv", "p_col": "p-value_missense", "exp_col": "expected_both_missense", "out": "df_top10_p_missense.csv"},
        {"file": "lof_final_result.csv", "p_col": "p-value_lof", "exp_col": "expected_both_lof", "out": "df_top10_p_lof.csv"},
        {"file": "lof_missense_combined_final_result.csv", "p_col": "p-value_lof_missense_combined", "exp_col": "expected_both_lof_missense_combined", "out": "df_top10_p_lof_missense.csv"},
       ]

    base_path = "/sci/labs/orzuk/mali.tsadok/UKB/singleton_1/results/final/lof/"  # change per folder if needed

    for cfg in configs:
        df = load_file(f"{base_path}{cfg['file']}")
        df = df[df[cfg['exp_col']] > 10]
        df = process_p_values(df, [cfg["p_col"]])
        df = df.sort_values(by=cfg["p_col"], ascending=True).reset_index(drop=True)
        df_top10 = df.nsmallest(10, cfg["p_col"]).reset_index(drop=True)
        df_top10.to_csv(f"/sci/labs/orzuk/mali.tsadok/UKB/singleton_1/results/final/{cfg['out']}", index=False)
