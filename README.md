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



# Cluster Execution Instructions

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

* The first phase uses data downloaded from DNA Nexus.
* Place this data in the **input folder** specified for the first phase in the CSV file.

## 3. Running the Pipeline

* For each phase, use the **corresponding shell script** from the CSV file.
* Each shell script calls a Python file, which is expected to be located at the path specified in the script.
* You can change the Python file location directly in the shell script if needed.
* Input and output folders for each phase can be modified in the CSV file under the **input** and **output** columns.


# Group-Based Permutation Pipeline

This pipeline performs **group-based permutation analysis** for two types of gene sets: `olida` and `paralogs`.

## Input Data
- Input files are located on the cluster.
- Make sure the data paths are updated in the scripts if needed.

## Pipeline Structure
The order of scripts is defined in the CSV `group_based_permutation_pipeline.csv`.  
Scripts are located in the corresponding folders in this repository:

- `olida/`
- `paralogs/`





