============
CAMP MAG QC
============


.. image:: https://readthedocs.org/projects/camp-mag_qc/badge/?version=latest
        :target: https://camp-mag_qc.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status

.. image:: https://img.shields.io/badge/version-0.1.0-brightgreen


Overview
--------

This module is designed to function as both a standalone MAG MAG QC pipeline as well as a component of the larger CAMP/CAP2 metagenome analysis pipeline. As such, it is both self-contained (ex. instructions included for the setup of a versioned environment, etc.), and seamlessly compatible with other CAMP modules (ex. ingests and spawns standardized input/output config files, etc.). 

The CAMP MAG quality-checking pipeline wraps several tools to assess the overall quality of the MAG binning process. The tools incorporated (CheckM, ``cmseq``, GTDB-Tk, ``dnadiff``) were inspired by Saheb Kashaf et al.'s MAG Snakemake workflow. More quality metrics will be added as they are validated on simulated and real gold-standard datasets. For example, adding QUAST's ability to identify contiguity and misassemblies would be useful.

Installation
------------

1. Clone repo from `Github <https://github.com/MetaSUB-CAMP/camp_mag_qc>`_. 

2. Set up the conda environment (contains, Snakemake) using ``configs/conda/mag_qc.yaml``. 

3. Make sure the installed pipeline works correctly. ``pytest`` only generates temporary outputs so no files should be created.
::
    cd camp_mag_qc
    conda env create -f configs/conda/mag_qc.yaml
    conda activate mag_qc
    pytest .tests/unit/

Using the Module
----------------

**Input**: ``/path/to/samples.csv`` provided by the user.
    * Note: The workflow expects MAG FastAs named with the format: ``{bin_num}.fa``, where the prefix ``bin_num`` is a unique identifier and the extension is ``fa`` (i.e.: standard MetaWRAP output). If this is not how your FastAs are formatted, either update the ``Snakefile`` or your MAG files.

**Output**: 1) An output config file summarizing 2) the module's outputs. 

- ``/path/to/work/dir/mag_qc/final_reports/{sample}.summary.csv`` that provides summary statistics of each MAG in the ``{sample}`` dataset:
    * Completeness contamination, as measured by CheckM's lineage-specific marker gene heuristic
    * Strain heterogeneity, as measured by ``cmseq``'s equation :math:`dN * 100 / dS`
    * Taxonomic classification, as estimated by GTDB-Tk's gene-based tree placement method
    * N50, MAG size (number of base-pairs), and GC content as adapted from a MetaWRAP script
    * MAG coverage by the closest reference genome, reference genome coverage by that MAG, average nucleotide identity between MAG and reference genome, as calculated by ``dnadiff``
- ``/path/to/work/dir/mag_qc/final_reports/samples.csv`` summarizing the locations of each of the summary reports

**Structure**:
::
    └── workflow
        ├── Snakefile
        ├── mag_qc.py
        ├── utils.py
        └── __init__.py
- ``workflow/mag_qc.py``: Click-based CLI that wraps the ``snakemake`` and unit test generation commands for clean management of parameters, resources, and environment variables.
- ``workflow/Snakefile``: The ``snakemake`` pipeline. 
- ``workflow/utils.py``: Sample ingestion and work directory setup functions, and other utility functions used in the pipeline and the CLI.

1. Make your own ``samples.csv`` based on the template in ``configs/samples.csv``. Sample test data can be found in ``test_data/``. 
    - ``ingest_samples`` in ``workflow/utils.py`` expects Illumina reads in FastQ (may be gzipped) form and de novo assembled contigs in FastA form
    - ``samples.csv`` requires either absolute paths or paths relative to the directory that the module is being run in

2. Update the relevant parameters in ``configs/parameters.yaml``.

3. Update the computational resources available to the pipeline in ``resources.yaml``. 

4. Download the database dependencies for GTDB-Tk and CheckM.
    * Note: ``GTDBTK_DATA_PATH`` has to be ``export``-ed every time you open a new instance of the terminal 
::

    export GTDBTK_DATA_PATH='/path/to/databases/metagenomics/GTDBTk_R202'
    wget https://data.gtdb.ecogenomic.org/releases/release202/202.0/auxillary_files/gtdbtk_r202_data.tar.gz -P ${GTDBTK_DATA_PATH}
    tar xvzf ${GTDBTK_DATA_PATH}/gtdbtk_r202_data.tar.gz -C ${GTDBTK_DATA_PATH} --strip 1

    export CHECK_DATA_PATH='/path/to/databases/metagenomics/CheckM_2015_01_16'
    wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz -P ${CHECK_DATA_PATH}
    tar xvzf ${CHECK_DATA_PATH}/checkm_data_2015_01_16.tar.gz -C ${CHECK_DATA_PATH} --strip 1

    # Tell CheckM where to find this data before running anything:
    checkm data setRoot ${CHECK_DATA_PATH}

5. To run CAMP on the command line, use the following, where ``/path/to/work/dir`` is replaced with the absolute path of your chosen working directory, and ``/path/to/samples.csv`` is replaced with your copy of ``samples.csv``. 
    - The default number of cores available to Snakemake is 1 which is enough for test data, but should probably be adjusted to 10+ for a real dataset.
    - Relative or absolute paths to the Snakefile and/or the working directory (if you're running elsewhere) are accepted!
::

    python /path/to/camp_mag_qc/workflow/mag_qc.py \
        (-c max_number_of_local_cpu_cores) \
        -d /path/to/work/dir \
        -s /path/to/samples.csv

* Note: This setup allows the main Snakefile to live outside of the work directory.

5. To run CAMP on a job submission cluster (for now, only Slurm is supported), use the following.
    - ``--slurm`` is an optional flag that submits all rules in the Snakemake pipeline as ``sbatch`` jobs. 
    - In Slurm mode, the ``-c`` flag refers to the maximum number of ``sbatch`` jobs submitted in parallel, **not** the pool of cores available to run the jobs. Each job will request the number of cores specified by threads in ``configs/resources/slurm.yaml``.
::

    sbatch -J jobname -o jobname.log << "EOF"
    #!/bin/bash
    python /path/to/camp_mag_qc/workflow/mag_qc.py --slurm \
        (-c max_number_of_parallel_jobs_submitted) \
        -d /path/to/work/dir \
        -s /path/to/samples.csv
    EOF

6. After checking over ``final_reports/`` and making sure you have everything you need, you can delete all intermediate files to save space. 
::

    python /path/to/camp_mag_qc/workflow/mag_qc.py cleanup \
        -d /path/to/work/dir \
        -s /path/to/samples.csv

7. If for some reason the module keeps failing, CAMP can print a script containing all of the remaining commands that can be run manually. 
::

    python /path/to/camp_mag_qc/workflow/mag_qc.py --dry_run \
        -d /path/to/work/dir \
        -s /path/to/samples.csv > cmds.txt
    python /path/to/camp_mag_qc/workflow/mag_qc.py commands cmds.txt

Extending the Module
Extending the Module
--------------------

We love to see it! This module was partially envisioned as a dependable, prepackaged sandbox for developers to test their shiny new tools in. 

These instructions are meant for developers who have made a tool and want to integrate or demo its functionality as part of the standard {{ cookiecutter.module_name }} workflow, or developers who want to integrate an existing tool. 

1. Write a module rule that wraps your tool and integrates its input and output into the pipeline. 
    - This is a great `Snakemake tutorial <https://bluegenes.github.io/hpc-snakemake-tips/>`_ for writing basic Snakemake rules.
    - If you're adding new tools from an existing YAML, use ``conda env update --file configs/conda/existing.yaml --prune``.
    - If you're using external scripts and resource files that i) cannot easily be integrated into either `utils.py` or `parameters.yaml`, and ii) are not as large as databases that would justify an externally stored download, add them to ``workflow/ext/`` or ``workflow/ext/scripts/`` and use ``rule external_rule`` as a template to wrap them. 
2. Update the ``make_config`` in ``workflow/Snakefile`` rule to check for your tool's output files. Update ``samples.csv`` to document its output if downstream modules/tools are meant to ingest it. 
    - If you plan to integrate multiple tools into the module that serve the same purpose but with different input or output requirements (ex. for alignment, Minimap2 for Nanopore reads vs. Bowtie2 for Illumina reads), you can toggle between these different 'streams' by setting the final files expected by ``make_config`` using the example function ``workflow_mode``.
    - Update the description of the ``samples.csv`` input fields in the CLI script ``workflow/{{ cookiecutter.module_slug }}.py``. 
3. If applicable, update the default conda config using ``conda env export > config/conda/{{ cookiecutter.module_slug }}.yaml`` with your tool and its dependencies. 
    - If there are dependency conflicts, make a new conda YAML under ``configs/conda`` and specify its usage in specific rules using the ``conda`` option (see ``first_rule`` for an example).
4. Add your tool's installation and running instructions to the module documentation and (if applicable) add the repo to your `Read the Docs account <https://readthedocs.org/>`_ + turn on the Read the Docs service hook.
5. Run the pipeline once through to make sure everything works using the test data in ``test_data/`` if appropriate, or your own appropriately-sized test data. Then, generate unit tests to ensure that others can sanity-check their installations.
    * Note: Python functions imported from ``utils.py`` into ``Snakefile`` should be debugged on the command-line first before being added to a rule because Snakemake doesn't port standard output/error well when using ``run:``.
::

    python /path/to/camp_mag_qc/workflow/mag_qc.py --unit_test \
        -d /path/to/work/dir \
        -s /path/to/samples.csv

6. Increment the version number of the modular pipeline.
::

    bump2version --allow-dirty --commit --tag major workflow/__init__.py \
                 --current-version A.C.E --new-version B.D.F

7. If you want your tool integrated into the main CAP2/CAMP pipeline, send a pull request and we'll have a look at it ASAP! 
    - Please make it clear what your tool intends to do by including a summary in the commit/pull request (ex. "Release X.Y.Z: Integration of tool A, which does B to C and outputs D").

.. ..

 <!--- 
 Bugs
 ----
 Put known ongoing problems here
 --->

Credits
-------

* This package was created with `Cookiecutter <https://github.com/cookiecutter/cookiecutter>`_ as a simplified version of the `project template <https://github.com/audreyr/cookiecutter-pypackage>`_.
* This module is heavily inspired by four Snakefiles from `MAG Snakemake workflow <https://github.com/Finn-Lab/MAG_Snakemake_wf>`_ (Saheb Kashaf et al. 2021).
* The MAG N50, size, and GC calculation rule was adapted from a script in `MetaWRAP <https://github.com/bxlab/metaWRAP>`_. 
* Free software: MIT 
* Documentation: https://mag-qc.readthedocs.io. 



