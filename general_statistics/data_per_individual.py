import os
from glob import glob
import pandas as pd
from collections import defaultdict
option_lof = [ 'frameshift_variant', 'frameshift_variant&splice_region_variant', 'frameshift_variant&start_lost',
               'frameshift_variant&stop_gained', 'frameshift_variant&stop_gained&splice_region_variant',
               'frameshift_variant&stop_lost', 'start_lost', 'stop_gained', 'stop_gained&splice_region_variant',
               'stop_lost', 'start_lost&conservative_inframe_deletion',
               'start_lost&conservative_inframe_insertion', 'stop_gained&disruptive_inframe_deletion',
               'stop_lost&conservative_inframe_deletion', 'start_lost&splice_region_variant',
               'stop_gained&disruptive_inframe_deletion&splice_region_variant', 'stop_lost&splice_region_variant',
               'bidirectional_gene_fusion', 'stop_gained&conservative_inframe_insertion',
               'stop_gained&disruptive_inframe_insertion', 'start_lost&disruptive_inframe_deletion',
               'stop_lost&disruptive_inframe_deletion', 'gene_fusion', 'rare_amino_acid_variant',
               'stop_gained&conservative_inframe_insertion&splice_region_variant',
               'stop_gained&disruptive_inframe_insertion&splice_region_variant',
               'frameshift_variant&synonymous_variant',
               'start_lost&conservative_inframe_deletion&splice_region_variant',
               'frameshift_variant&stop_lost&splice_region_variant',
               'start_lost&disruptive_inframe_insertion',
               'frameshift_variant&missense_variant&splice_region_variant',
               'frameshift_variant&start_lost&splice_region_variant',
               'frameshift_variant&missense_variant',
               'stop_lost&conservative_inframe_deletion&splice_region_variant',
               'stop_gained&conservative_inframe_deletion',
               'start_lost&disruptive_inframe_deletion&splice_region_variant',
               'stop_lost&disruptive_inframe_deletion&splice_region_variant',
               'start_lost&conservative_inframe_insertion&splice_region_variant',
               'stop_lost&disruptive_inframe_insertion' ]

option_missense = ["missense_variant", "missense_variant&splice_region_variant",
                   "missense_variant&disruptive_inframe_insertion",
                   "missense_variant&conservative_inframe_insertion"]  # adjust if needed

def process_all_subfolders(root_dir):
    all_files = glob(os.path.join(root_dir, '**/*.csv'), recursive=True)

    results = defaultdict(lambda: {
        "all_genes": set(),
        "lof_genes": set(),
        "missense_genes": set(),
        "all_genes_homo": set(),
        "lof_genes_homo": set(),
        "missense_genes_homo": set()
    })

    for file_path in all_files:
        df = pd.read_csv(file_path)

        print (
        
        # all genes
        all_counts = df.groupby("Individual")["gene"].unique()
        for ind, genes in all_counts.items():
            results[ind]["all_genes"].update(genes)

        # LOF
        lof_counts = df[df["mutation_type"].isin(option_lof)].groupby("Individual")["gene"].unique()
        for ind, genes in lof_counts.items():
            results[ind]["lof_genes"].update(genes)

        # Missense
        missense_counts = df[df["mutation_type"].isin(option_missense)].groupby("Individual")["gene"].unique()
        for ind, genes in missense_counts.items():
            results[ind]["missense_genes"].update(genes)

        # Homozygous
        df_homo = df[df["count_Mut"] == 2]

        all_homo_counts = df_homo.groupby("Individual")["gene"].unique()
        for ind, genes in all_homo_counts.items():
            results[ind]["all_genes_homo"].update(genes)

        lof_homo_counts = df_homo[df_homo["mutation_type"].isin(option_lof)].groupby("Individual")["gene"].unique()
        for ind, genes in lof_homo_counts.items():
            results[ind]["lof_genes_homo"].update(genes)

        missense_homo_counts = df_homo[df_homo["mutation_type"].isin(option_missense)].groupby("Individual")["gene"].unique()
        for ind, genes in missense_homo_counts.items():
            results[ind]["missense_genes_homo"].update(genes)
        break

    # Convert results to DataFrame
    final_data = []
    for ind, vals in results.items():
        final_data.append({
            "Individual": ind,
            "all_genes": len(vals["all_genes"]),
            "lof_genes": len(vals["lof_genes"]),
            "missense_genes": len(vals["missense_genes"]),
            "all_genes_homo": len(vals["all_genes_homo"]),
            "lof_genes_homo": len(vals["lof_genes_homo"]),
            "missense_genes_homo": len(vals["missense_genes_homo"])
        })

    return pd.DataFrame(final_data)


if __name__ == "__main__":
    result =  process_all_subfolders("/sci/labs/orzuk/mali.tsadok/UKB/singleton_1/filter_data/lof")
    result.to_csv("/sci/labs/orzuk/mali.tsadok/UKB/data_per_individual_new_model_1_percent.csv" , index= False)