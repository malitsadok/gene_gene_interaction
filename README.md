# Pair gens analysis 

## DNA Nexus Part
This part of the project is handled by the folder **annotation_applet**.  
It contains the necessary files and structure to run the annotation workflow in DNA Nexus.

## Requirements
- DNA Nexus user account.
- Create an applet on DNA Nexus with the **same structure** as the `annotation_applet` folder.
- Each DNA Nexus applet automatically runs the `code.sh` script when executed.

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

* `--destination` specifies the output folder.
* The batch TSV file contains all VCF files for chromosome 16.

The batch TSV contains the paths for all files in chromosome 16.



## Cluster Part

This part of the project is designed to run on the **university cluster**.
It uses the folder **pipeline_scripts_singleton** and executes the pipeline in multiple phases.

### Requirements

* Access to the university cluster.
* The `pipeline_scripts_singleton` folder containing all necessary scripts.
* A CSV file specifying the phases, which includes:

  * The shell script to run for each phase.
  * The input folder for that phase.
  * The output folder for that phase.

### Workflow

1. **Phase Execution:**
   Each phase must be run **one by one** according to the order specified in the CSV file.

2. **Input Data:**

   * The first phase uses data downloaded from DNA Nexus.
   * This data should be placed in the input folder of the first phase.

3. **Running Phases:**

   * For each phase, use the corresponding shell script from the CSV.
   * Each shell script calls a Python file; the script expects the Python file to be located in the path it is looking for.
   * For the singleton pipeline, both the shell and Python files are in the `singleton_pipeline` folder.
   * Make sure the input and output folders are correctly specified for each run.

### Running the Pipeline on All-Set Data

* The same pipeline can be run on the **all-set version of the data**.
* For this version, there are some differences in the scripts and in the folder structure for storing the data.
* The new pipeline is defined in the CSV file `pipeline_scripts_all_data`  and the shell and Python files are in the `all_set_pipeline` folder.



