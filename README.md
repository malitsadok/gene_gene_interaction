# Pair genes pipeline

## DNA Nexus Part
This part of the project is handled by the folder **annotation_applet**.  
It contains the necessary files and structure to run the annotation workflow in DNA Nexus.

## Requirements
- DNA Nexus user account.
- Create an applet on DNA Nexus with the **same structure** as the `annotation_applet` folder.
- Each DNA Nexus applet automatically runs the `code.sh` script when executed.

example of how running the applet  chromosome 21 block number 4 :
dx run annotation_applet -i vcf_in='Bulk/Exome sequences/Population level exome OQFE variants, pVCF format - final release/ukb23157_c21_b4_v1.vcf.gz' --destination "/results/"

## code.sh Details
The `code.sh` script runs a public Docker image:

```bash
docker run -v $PWD/analysis:/analysis malitsadok/image2:1
````

* The Docker image is public; anyone can access it.
* Processes one VCF file at a time, adds annotations and filtering, and returns a CSV file.

## Input Data

* The DNA Nexus project contains ~1,000 VCF files.
* Each run processes a single VCF file, adds annotations, applies filtering, and produces a CSV output.

## Running in Parallel / Batch Mode

To run multiple files in parallel, use **batch mode**.
Example for chromosome 16:

```bash
dx run annotation_applet --batch-tsv thesa/UKbiobank/batch_files/chr16.0000.tsv --destination="Gene Pair Interactions:/results/chr16"
```
 explanation how to create batch file you can find here:
 https://documentation.dnanexus.com/user/running-apps-and-workflows/running-batch-jobs



## Cluster Execution Instructions

This project includes two pipelines: **singleton_pipeline** and **all_data_pipeline**.  

* **all_data_pipeline**: generates gene-pair combinations considering all variants (both singleton and non-singleton).  
* **singleton_pipeline**: generates gene-pair combinations where at least one of the genes is restricted to singleton variants.  

Each pipeline contains:

* A CSV file (`scripts_singleton_pipeline.csv` or `scripts_all_data_pipeline.csv`) specifying the phases, which includes:
  * The shell script to run for each phase.
  * The input folder for that phase.
  * The output folder for that phase.

## 1. Phase Execution

Each phase must be run **one by one**, following the order specified in the CSV file.

## 2. Input Data

* The first phase uses data downloaded from DNA Nexus. and can be fount in
?????????????????????????????
* Place this data in the **input folder** specified for the first phase in the CSV file.

## 3. Running the Pipeline

* For each phase, use the **corresponding shell script** from the CSV file.
* Each shell script calls a Python file, which is expected to be located at the path specified in the script.
* You can change the Python file location directly in the shell script if needed.
* Input and output folders for each phase can be modified in the CSV file under the **input** and **output** columns.


Sure — here is a **clean, organized Markdown version** you can paste directly into your `README.md`. I’ve kept the content faithful but improved structure, clarity, and grammar.

---

# Group-Based Permutation Pipeline

This pipeline performs **group-based permutation analysis** for two types of gene sets:

* **Olida**
* **Paralogs**



### Shared Input Files

The `all_data_pipeline` generates the following ZIP archive, which serves as the **input** for this pipeline:

```
all_data_results.zip
```

This archive contains three CSV files:

* `lof_final_result.csv`
* `missense_final_result.csv`
* `lof_missense_combined_final_result.csv`

These files are required for all downstream analyses in both the **Olida** and **Paralogs** pipelines.


## Olida Pipeline

### Step 1: Preparation

Run the Python script:

```
olida_preperations.py
```

#### Input Files

* `olida.json`
* `GeneCombination_olida.csv`

#### Output Files

This step generates three CSV files, split by score category:

* `olida_pairs_filtered_by_score_0_original.csv`
* `olida_pairs_filtered_by_score_1_original.csv`
* `olida_pairs_filtered_by_score_2_3_original.csv`



### Step 2: Run Score-Specific Analysis (Example: Score = 0)

The following steps must be run **separately for each score category**.

Below, we demonstrate the workflow using **score = 0** as an example.

> ⚠️ **Important:**
> For each score category (e.g., score 0, 1, or 2–3), you must manually update the input file paths and file names in the scripts below so that they correspond to the correct score.

---

### Step 2.1: Add Expected and Observed Counts

Run:

```
add_expected_for_olida.py
```

#### Inputs

* `olida_pairs_filtered_by_score_0_original.csv`
* The three CSV files from `all_data_results.zip`

#### Output

* Adds **expected** and **observed** counts for score 0 gene pairs

---

### Step 2.2: Generate Permutations

Run:

```
create_permutation_olida_model.py
```

#### Description

* Generates **10,000 permutations** for the Olida model
* Produces permutation files for score 0

---

### Step 2.3: Add Expected Values to Permutations

Run:

```
add_expected_for_olida_permutations.py
```

#Perfect — here is the **corrected version** that makes it clear this step is **run once** and **aggregates all score categories**, while keeping the same style as the rest of your README.

---

Got it — you don’t want this to be labeled as part of **Step 2** at all.
Here is the **correct fix**: renumber it as a **separate step** and clearly state it runs **once for all scores**.

You can replace the section with this:

---

## Step 3: Null Distribution and P-Value Calculation

After completing all permutations **for all score categories**, run:

```
olida_distribution.py
```

#### Description

* Constructs the **null distribution** using permutation results from **all score categories**
* Extracts **p-values** for the observed statistics

> ⚠️ **Important:**
> This step is executed **once**, after all score-specific analyses are complete.
> Ensure that folder paths are correctly configured so the script can access:
>
> * Original (observed) files for all scores
> * Permutation files for all scores


