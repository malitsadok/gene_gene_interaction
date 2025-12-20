import os
import pandas as pd
import sys


MISSENSE_OPTIONS = [
    "missense_variant",
    "missense_variant&splice_region_variant",
    "missense_variant&disruptive_inframe_insertion",
    "missense_variant&conservative_inframe_insertion"
]


def is_missense(mutation_type):
    if pd.isna(mutation_type):
        return False
    return any(opt in mutation_type for opt in MISSENSE_OPTIONS)


def count_appearance_for_each_snv(df):
    # Count the number of individuals who have each SNV
    snv_count_df = pd.DataFrame({
        "count": df.groupby(
            ['SNV', 'count_Mut', 'gene', 'mutation_type']
        ).count_Mut.count()
    })

    snv_count_df = (
        snv_count_df
        .pivot_table(
            index=['SNV', 'gene', 'mutation_type'],
            columns='count_Mut',
            values='count',
            fill_value=0
        )
        .reset_index()
    )

    snv_count_df.rename(columns={1: 'countHetero', 2: 'countHomo'}, inplace=True)
    return snv_count_df


def process_file(file_path):
    print("Processing:", file_path)
    df = pd.read_csv(file_path)

    snv_count_df = count_appearance_for_each_snv(df)

    columns = ['Individual', 'gene', 'count_Mut']

    # Identify singleton variants
    singleton_vars = snv_count_df[
        (snv_count_df.countHetero == 1) &
        (snv_count_df.countHomo == 0)
    ]

    singleton_vars_missense = singleton_vars[
        singleton_vars['mutation_type'].apply(is_missense)
    ]

    singleton_vars_lof = singleton_vars[
        ~singleton_vars['mutation_type'].apply(is_missense)
    ]

    is_missense_series = df['mutation_type'].apply(is_missense)

    categories = {
        'singleton_lof': (
            df.SNV.isin(singleton_vars_lof.SNV),
            0, 1, 0, 0
        ),
        'singleton_missense': (
            df.SNV.isin(singleton_vars_missense.SNV),
            0, 0, 0, 1
        ),
        'missense': (
            is_missense_series,
            0, 0, 1, 0
        ),
        'lof': (
            ~is_missense_series,
            1, 0, 0, 0
        ),
    }

    dfs = []

    for category, (condition, lof, lof_singleton, missense, missense_singleton) in categories.items():
        filtered_df = df[condition]

        filtered_no_dup_df = (
            filtered_df
            .drop_duplicates(subset=['gene', 'Individual'])[columns]
            .copy()
        )

        filtered_no_dup_df["lof"] = lof
        filtered_no_dup_df["lof_singleton"] = lof_singleton
        filtered_no_dup_df["missense"] = missense
        filtered_no_dup_df["missense_singleton"] = missense_singleton

        dfs.append(filtered_no_dup_df)

    return pd.concat(dfs, ignore_index=True)


def process(input_directory, output_directory):
    files = [
        f for f in os.listdir(input_directory)
        if os.path.isfile(os.path.join(input_directory, f))
    ]

    print("Files:", files)

    os.makedirs(output_directory, exist_ok=True)

    for file_name in files:
        full_path = os.path.join(input_directory, file_name)

        result_df = process_file(full_path)
        result_df = result_df.drop("count_Mut", axis=1)

        grouped_df = result_df.groupby(
            ['Individual', 'gene'],
            as_index=False
        ).max()

        grouped_df.to_csv(
            os.path.join(output_directory, file_name),
            index=False
        )


if __name__ == "__main__":
    chr_name = sys.argv[1]
    input_directory = sys.argv[2]
    output_directory = sys.argv[3]

    print("Chromosome:", chr_name)

    process(
        input_directory + chr_name,
        output_directory + chr_name
    )
