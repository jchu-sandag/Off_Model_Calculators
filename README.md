# Off_Model_Calculators

Off-model calculators to capture VMT and GHG reduction from vanpool and carshare strategies not addressed by ABM3

If running the process for the first time, run all steps to create the conda environment. If repeating the process, simply run step 2 and 5

1. To create the conda environment for OMC, run
conda create -n env_omc

______________________________________________________________________________________________
2. To activate conda environment, run
conda activate env_omc

______________________________________________________________________________________________
3. To install packages to conda environment, run
conda install --file=env_packages.txt

______________________________________________________________________________________________
4. To install "openmatrix" package to conda environment, run
conda install conda-forge::openmatrix

______________________________________________________________________________________________
5. To run the calculators in batch format, open cmd in the directory corresponding
to ".bat" file and call "omc_batch-abm3.bat"

______________________________________________________________________________________________
** Note: process developed based on conda version 24.3.0. The command for other versions might be different


