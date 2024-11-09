# MAG QC

[![Documentation Status](https://img.shields.io/badge/docs-passing-brightgreen.svg)](https://camp-documentation.readthedocs.io/en/latest/magqc.html) ![Version](https://img.shields.io/badge/version-0.13.0-brightgreen)

<!-- [![Documentation Status](https://img.shields.io/readthedocs/camp-mag_qc)](https://camp-documentation.readthedocs.io/en/latest/mag_qc.html) -->

## Overview

This module is designed to function as both a standalone MAG QC pipeline as well as a component of the larger CAMP metagenome analysis pipeline. As such, it is both self-contained (ex. instructions included for the setup of a versioned environment, etc.), and seamlessly compatible with other CAMP modules (ex. ingests and spawns standardized input/output config files, etc.). 

The CAMP MAG quality-checking pipeline wraps several tools (CheckM, gunc, GTDB-Tk, DNADiff, QUAST) to assess the overall quality of the MAG binning process. More quality metrics will be added as they are validated on simulated and real gold-standard datasets. 

## Installation

> [!TIP]
> All databases used in CAMP modules will also be available for download on Zenodo (link TBD).

1. Clone repo from [Github](<https://github.com/MetaSUB-CAMP/camp_mag-qc>).
```Bash
git clone https://github.com/MetaSUB-CAMP/camp_mag-qc
```

2. Set up the conda environment using `configs/conda/mag-qc.yaml`. 

If you don't already have `conda` handy, we recommend installing `miniforge`, which is a minimal conda installer that, by default, installs packages from open-source community-driven channels such as `conda-forge`.
```Bash
# If you don't already have conda on your system...
# wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh

# Create and activate conda environment 
cd camp_mag-qc
conda env create -f configs/conda/mag_qc.yaml
conda activate mag-qc
```

3. Download the database dependencies for GTDB-Tk.
    - Note: `GTDBTK_DATA_PATH` has to be `export`-ed every time you open a new instance of the terminal 
```Bash
export GTDBTK_DATA_PATH='/path/to/databases/metagenomics/GTDBTk_R202'
wget https://data.gtdb.ecogenomic.org/releases/release202/202.0/auxillary_files/gtdbtk_r202_data.tar.gz -P ${GTDBTK_DATA_PATH}
tar xzf ${GTDBTK_DATA_PATH}/gtdbtk_r202_data.tar.gz -C ${GTDBTK_DATA_PATH} --strip 1
```

4. Download the database dependencies for CheckM2. The easiest way to do this is to install the CheckM2 environment using `--dry_run` (see below for explanation) and then activating the CheckM2 conda environment to use its database download command.
```Bash
python /path/to/camp_mag-qc/workflow/mag_qc.py --dry_run \
    -d /home/lam4003/bin/camp_mag-qc/test_out \
    -s /home/lam4003/bin/camp_mag-qc/test_data/samples.csv

# In the directory /path/to/camp_mag-qc/conda_envs/, find the environment ID that corresponds to CheckM2

# Activate the conda environment that corresponds to CheckM2
conda activate /path/to/camp_mag-qc/conda_envs/checkm2_env_id

checkm2 database --download --path /path/to/databases
```

5. Download the database dependencies for CheckM1. 
```
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz -P /path/to/databases
tar xzf checkm_data_2015_01_16.tar.gz
```

5. Download the database dependencies for gunc. 
```Bash
gunc download_db /path/to/databases
```

6. Update the parameters `ext` and the paths to the databases for GTDB-Tk, CheckM1/2, gunc, and DIAMOND in `test_data/parameters.yaml`.

7. Make sure the installed pipeline works correctly. With 40 threads and a maximum of 250 GB allocated, the test dataset should finish in approximately 43 minutes. `pplacer`, a component of GTDB-Tk, requires a large amount of RAM to complete.
```Bash
# Run tests on the included sample dataset
python /path/to/camp_mag-qc/workflow/mag_qc.py test
```

## Using the Module

**Input**: `/path/to/samples.csv` provided by the user.
    - Note that the `samples.csv` generated by the CAMP binning module needs to be manually pruned so that there is only a single column called 'mag_dir' indicating the absolute path to the directory containing the MAG FastA files.
    - The workflow expects MAG FastAs named with the format: `bin.{bin_num}.fa`, where the prefix `bin_num` is a unique identifier and the extension is `fa` (i.e.: standard MetaWRAP output). If this is not how your FastAs are formatted, either update `ingest_samples` in `workflow/utils.py` or your MAG files.

**Output**: 1) An output config file summarizing 2) the module's outputs. See `test_data/test_out.tar.gz` for a sample output work directory.

- `/path/to/work/dir/mag_qc/final_reports/{sample}.summary.csv` that provides summary statistics of each MAG in the `{sample}` dataset:
    * Completeness contamination, N50, MAG size (number of base-pairs), and GC content as measured by CheckM's lineage-specific marker gene heuristic
    * Taxonomic classification, as estimated by GTDB-Tk's gene-based tree placement method
    * MAG coverage by the closest reference genome, reference genome coverage by that MAG, average nucleotide identity between MAG and reference genome, as calculated by `dnadiff`
    * Reference genome-based metrics such as NG50, NA50, proportion of misassembled contigs and sequence data in misassemblies, etc.

### Module Structure
```
└── workflow
    ├── Snakefile
    ├── mag_qc.py
    ├── utils.py
    ├── __init__.py
    └── ext/
        └── scripts/
```
- `workflow/mag_qc.py`: Click-based CLI that wraps the `snakemake` and other commands for clean management of parameters, resources, and environment variables.
- `workflow/Snakefile`: The `snakemake` pipeline. 
- `workflow/utils.py`: Sample ingestion and work directory setup functions, and other utility functions used in the pipeline and the CLI.
- `ext/`: External programs, scripts, and small auxiliary files that are not conda-compatible but used in the workflow.

### Running the Workflow

1. Make your own `samples.csv` based on the template in `configs/samples.csv`. Sample test data can be found in `test_data/`. 
    - `samples.csv` requires either absolute paths or paths relative to the directory that the module is being run in

2. Update the relevant parameters in `configs/parameters.yaml`.

3. Update the computational resources available to the pipeline in `configs/resources.yaml`. 

#### Command Line Deployment

To run CAMP on the command line, use the following, where `/path/to/work/dir` is replaced with the absolute path of your chosen working directory, and `/path/to/samples.csv` is replaced with your copy of `samples.csv`. 
    - The default number of cores available to Snakemake is 1 which is enough for test data, but should probably be adjusted to 10+ for a real dataset.
    - Relative or absolute paths to the Snakefile and/or the working directory (if you're running elsewhere) are accepted!
    - The parameters and resource config YAMLs can also be customized.
```Bash
python /path/to/camp_mag-qc/workflow/mag_qc.py \
    (-c number_of_cores_allocated) \
    (-p /path/to/parameters.yaml) \
    (-r /path/to/resources.yaml) \
    -d /path/to/work/dir \
    -s /path/to/samples.csv
```

#### Slurm Cluster Deployment

To run CAMP on a job submission cluster (for now, only Slurm is supported), use the following.
    - `--slurm` is an optional flag that submits all rules in the Snakemake pipeline as `sbatch` jobs. 
    - In Slurm mode, the `-c` flag refers to the maximum number of `sbatch` jobs submitted in parallel, **not** the pool of cores available to run the jobs. Each job will request the number of cores specified by threads in `configs/resources/slurm.yaml`.
```Bash
sbatch -J jobname -o jobname.log << "EOF"
#!/bin/bash
python /path/to/camp_mag-qc/workflow/mag_qc.py --slurm \
    (-c max_number_of_parallel_jobs_submitted) \
    (-p /path/to/parameters.yaml) \
    (-r /path/to/resources.yaml) \
    -d /path/to/work/dir \
    -s /path/to/samples.csv
EOF
```

#### Finishing Up

1. To plot all of the metrics extracted above, set up the dataviz environment and follow the instructions in the Jupyter notebook:
```Bash
conda env create -f configs/conda/dataviz.yaml
conda activate dataviz
jupyter notebook &
```

2. After checking over `final_reports/` and making sure you have everything you need, you can delete all intermediate files to save space. 
```Bash
python3 /path/to/camp_mag-qc/workflow/mag_qc.py cleanup \
    -d /path/to/work/dir \
    -s /path/to/samples.csv
```

3. If for some reason the module keeps failing, CAMP can print a script containing all of the remaining commands that can be run manually. 
```Bash
python3 /path/to/camp_mag-qc/workflow/mag_qc.py --dry_run \
    -d /path/to/work/dir \
    -s /path/to/samples.csv
```

## Credits

- This package was created with [Cookiecutter](https://github.com/cookiecutter/cookiecutter>) as a simplified version of the [project template](https://github.com/audreyr/cookiecutter-pypackage>).
- Free software: MIT

