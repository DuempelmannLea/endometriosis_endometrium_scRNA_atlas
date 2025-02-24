# Snakemake-seq workflow
This is a pipeline for analyzing single-cell RNA sequencing data. Input files are either fastq files or count matrices of RNA expression data.

## Setting up the snakemake workflow
Specific preparations required to run the pipeline:

1. Download the reference files (genome annotation file and the reference transcriptome) into the directory specified in the config file.

2. Edit the config file appropriately in the Snakemake working directory. The options are described in the config file.

3. Edit the sample file, it requires a comma-separated .csv file containing the following columns:

   - **names**: The names (without path) of the gzipped fastq files. One name should correspond to one original biological sample (excluding lane and read information), e.g., specify M\_15\_S15 for file M\_15\_S15\_L002\_R2\_001.fastq.gz.&#x20;
   - **type**: Defines the fragment ends from where the sequencing was performed. Either PE (paired end) or SE (single end).
   - **celline**: The name of the cell line or cell type that was sequenced (e.g., PBMC, or NIH-3T3).
   - **treatment**: The treatment that this sample received. This is typically used as a grouping variable to create contrasts downstream.
   - **batch**: The batch for this file. This column is used to specify batch effects in downstream analyses.

   
1. Create a conda environment to run the pipeline:

   ```
   conda create --name ENV python=3.8
   conda activate ENV
   conda install mamba -n base -c conda-forge
   conda install pandas
   pip install snakemake==6.6
   ```

2. Run the pipeline, preferably by setting up a profile config.yaml and running it using:

   ```
   snakemake --profile <profile_repository>
   ```


## Running snakemake

A number of features ensure that the workflow is automatically and dynamically created based on the config file:

- The pipelines that need to be run are automatically determined as the targets from rule `all` in the main pipeline are dynamically created based on the config parameters.
- All files and the metadata contained in their names (such as read/lane information) will be automatically detected and used within the workflow. (The location of the files is defined via `fastq_dir` in the config file.)

Snakemake can be run from the server or in the cloud. The following flags can be added to a run:

- `-F` forces the workflow to rerun.
- `--dry-run` performs a dry run to test dependencies between rules.
- `-nt` ensures temporary files will not be deleted, which is useful for inspection and debugging.
- `-k` allows independent jobs to continue running after other jobs fail.
- `--dag | dot -Tpdf > dag.pdf` or `--report` generates a directed acyclic graph showing rule dependencies.
- `-d` sets the root directory for logging files.

### Running using a profile

Use a configuration profile where default options for the Snakemake CLI are stored. A template profile can be found in `/templates/profile`. Snakemake can then be invoked using:

```
snakemake --profile profile
```

Alternatively, an absolute or relative path to the private profile directory can be given.

## Post-run

### Output

The exact output of a run depends on the configured components. The files are located:

- **In the user-defined output folder**:

  - Intermediate output generally consists of binary `.rds` files.
  - Logs contain a transcript from the command line as `.log` files and a copy of the R environment as `.RData` files.
  - Reports are `.html` files.
  - Plots are intermediate plots used to generate reports.

- **In the folder from which the pipeline was run**:

  - A log file if output was redirected as described above.

### Cleanup

After a successful run, remove the conda environments using:

```
snakemake --cleanup-conda
```

## Test data

FASTQ files for testing and debugging can be downloaded from 10X by running:

```
wget http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar
```

The data consist of two sample PBMCs (pbmc\_1k\_v3\_S1\_L001 and pbmc\_1k\_v3\_S1\_L001).

## Debugging

If the Alevin quantification throws an error, check the number of cells found by Alevin. Try running Alevin with the `--expectCells X` flag, where X is the expected number of cells. This skips knee calculation and runs Alevin in a CellRanger-like mode. For errors occurring in R scripts, debugging should be done directly using Rscript. A `snakemake` object can be loaded for convenience, containing all relevant parameters (input and output paths, config parameters, etc.). The parsed objects are located in:

```
<output>/logs/<rule>/<rulename>.rds
```

