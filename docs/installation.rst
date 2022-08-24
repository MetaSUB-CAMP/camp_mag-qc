.. highlight:: shell

============
Installation
============


Stable release
--------------

1. Clone repo from `github <https://github.com/lauren-mak/camp_mag_qc>_`. 

2. Set up the conda environment (contains, Snakemake) using ``configs/conda/camp_mag_qc.yaml``. 

3. Make sure the installed pipeline works correctly. ``pytest`` only generates temporary outputs so no files should be created.
::
    cd camp_mag_qc
    conda env create -f configs/conda/camp_mag_qc.yaml
    conda activate camp_mag_qc
    pytest .tests/unit/

